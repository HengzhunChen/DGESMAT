function HamDG = Setup(HamDG, esdfParam)
% HAMILTONIANDG/SETUP initializes HamiltonianDG object HamDG with data 
%    from ESDFInputParam object esdfParam.
% 
%    See also HamiltonianDG, Domain, Atom, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% **********************************************************************
% basic parameters
% **********************************************************************

HamDG.domain            = esdfParam.basic.domain;
HamDG.atomList          = esdfParam.basic.atomList;
HamDG.pseudoType        = esdfParam.basic.pseudoType;
HamDG.XCType            = esdfParam.basic.XCType;
HamDG.numExtraState     = esdfParam.basic.numExtraState;
HamDG.ecutWavefun       = esdfParam.basic.ecutWavefunction;
HamDG.numElem           = esdfParam.DG.numElem;
HamDG.penaltyAlpha      = esdfParam.DG.penaltyAlpha;

% NOTE: only consider numSpin == 2 in the DG calculation
HamDG.numSpin = 2;

for d = 1 : dimDef()
    if mod( HamDG.domain.numGrid(d), HamDG.numElem(d) ) ~= 0
        error('The number of global wfc grid points is not divisible by the number of elements');
    end
    if mod( HamDG.domain.numGridFine(d), HamDG.numElem(d) ) ~= 0
        error('The number of global rho grid points is not divisible by the number of elements');
    end
end

% fft over global domain
HamDG.fft = Fourier(esdfParam.basic.domain);

% DG Hamiltonian matrix
numElemTotal = prod(HamDG.numElem);
HamDG.HMat = cell(numElemTotal, numElemTotal);
HamDG.sizeHMat = sum(esdfParam.DG.numALBElem);



% ***********************************************************************
%  Partition of elements
% ***********************************************************************

numElem = HamDG.numElem;
numElemTotal = prod(numElem);

HamDG.grid.numLGLGridElem    = esdfParam.DG.numGridLGL;
HamDG.grid.numUniformGridElem = HamDG.domain.numGrid ./ numElem;
HamDG.grid.numUniformGridElemFine = HamDG.domain.numGridFine ./ numElem;

% Setup the element domains
HamDG.domainElem = cell(numElemTotal, 1);
for elemIdx = 1 : numElemTotal
    [i, j, k] = ElemIdxToKey(elemIdx, numElem);
    key = [i-1, j-1, k-1];
    dm = Domain();
    dm.length = HamDG.domain.length ./ numElem;
    dm.numGrid = HamDG.grid.numUniformGridElem;
    dm.numGridFine = HamDG.grid.numUniformGridElemFine;
    dm.posStart = dm.length .* key;
    HamDG.domainElem{elemIdx} = dm;
end


% Partition by element
HamDG.density     = cell(numElemTotal, 1);
HamDG.densityLGL  = cell(numElemTotal, 1);
HamDG.vext        = cell(numElemTotal, 1);
HamDG.vhart       = cell(numElemTotal, 1);
HamDG.vxc         = cell(numElemTotal, 1);
HamDG.epsxc       = cell(numElemTotal, 1);
HamDG.vtot        = cell(numElemTotal, 1);
HamDG.vtotLGL     = cell(numElemTotal, 1);

% Initialize quantities
for elemIdx = 1 : numElemTotal
    numGridElemFineTotal = prod(HamDG.grid.numUniformGridElemFine);
    numLGLGridElemTotal  = prod(HamDG.grid.numLGLGridElem);
    
    HamDG.density{elemIdx}     = zeros(numGridElemFineTotal, 1);
    HamDG.densityLGL{elemIdx}  = zeros(numLGLGridElemTotal, 1);
    HamDG.vext{elemIdx}        = zeros(numGridElemFineTotal, 1);
    HamDG.vhart{elemIdx}       = zeros(numGridElemFineTotal, 1);
    HamDG.vxc{elemIdx}         = zeros(numGridElemFineTotal, 1);
    HamDG.epsxc{elemIdx}       = zeros(numGridElemFineTotal, 1);
    HamDG.vtot{elemIdx}        = zeros(numGridElemFineTotal, 1);
    HamDG.vtotLGL{elemIdx}     = zeros(numLGLGridElemTotal, 1);
end

HamDG.gradDensity = cell(dimDef(), 1);
for d = 1 : dimDef()
    HamDG.gradDensity{d} = cell(numElemTotal, 1);
    for elemIdx = 1 : numElemTotal
        HamDG.gradDensity{d}{elemIdx} = zeros(numGridElemFineTotal, 1);
    end
end



% **********************************************************************
% Generate the grids
% **********************************************************************

% When dual grid is used, the fine grid is used for quantities such as
% density and potential. The coarse grid is used for wavefunctions (basis
% functions)

HamDG.grid.uniformGrid = UniformMesh(HamDG.domain);
HamDG.grid.uniformGridFine = UniformMeshFine(HamDG.domain);

HamDG.grid.uniformGridElem = cell(numElemTotal, 1);
HamDG.grid.uniformGridElemFine = cell(numElemTotal, 1);
HamDG.grid.LGLGridElem = cell(numElemTotal, 1);

for elemIdx = 1 : numElemTotal
    HamDG.grid.uniformGridElem{elemIdx} = ...
        UniformMesh( HamDG.domainElem{elemIdx} );
    HamDG.grid.uniformGridElemFine{elemIdx} = ...
        UniformMeshFine( HamDG.domainElem{elemIdx} );
    HamDG.grid.LGLGridElem{elemIdx} = ...
        LGLMesh( HamDG.domainElem{elemIdx}, HamDG.grid.numLGLGridElem );
end


% ----------------------------------------------------------------------
% Generate the differentiation matrix on the LGL grid
% NOTE: This assumes uniform mesh used for each element
% ----------------------------------------------------------------------

for d = 1 :  dimDef()
    [~, ~, ~, DMat] = GenerateLGL( HamDG.grid.numLGLGridElem(d) );    
    % Scale the differentiation matrix
    DMat = DMat .* ( 2 ./ (HamDG.domain.length(d) ./ numElem(d)) ); 
    HamDG.grid.DMat{d} = DMat;
end


% --------------------------------------------------------------------- 
% Transfer matrix from LGL grid to uniform grid 
% ---------------------------------------------------------------------

% Generate the transfer matrix from LGL grid to uniform grid on each
% element. The Lagrange polynomials involved in the transfer matrix is
% computed using the Barycentric method. For more information see
%
% [J.P. Berrut and L.N. Trefethen, Barycentric Lagrange Interpolation,
% SIAM Rev. 2004]
%
% NOTE: This assumes uniform mesh used for each element.

domainLength = HamDG.domainElem{1}.length;
oldGrid = HamDG.grid.LGLGridElem{1};

newGrid = HamDG.grid.uniformGridElem{1};
HamDG.grid.LGLToUniformMat = GenerateLagrangeInterpMat(...
        newGrid, oldGrid, domainLength);

newGrid = HamDG.grid.uniformGridElemFine{1};
HamDG.grid.LGLToUniformMatFine = GenerateLagrangeInterpMat(...
        newGrid, oldGrid, domainLength);


% ----------------------------------------------------------------------
% Compute the LGL weights at 1D, 2D, 3D
% ----------------------------------------------------------------------

dmElemlength = HamDG.domainElem{1}.length;
numGrid = HamDG.grid.numLGLGridElem;

% compute the integration weights
% 1D
for d = 1 : dimDef()
    [~, LGLWeight1D, ~, ~] = GenerateLGL(numGrid(d));
    HamDG.grid.LGLWeight1D{d} = 0.5 * dmElemlength(d) .* LGLWeight1D;
end

LGLWeight1D = HamDG.grid.LGLWeight1D;

% 2D: faces labeled by normal vectors, i.e.
% yz face : 1
% xz face : 2
% xy face : 3

% yz face
HamDG.grid.LGLWeight2D{1} = (LGLWeight1D{2})' .* LGLWeight1D{3}; 
% xz face
HamDG.grid.LGLWeight2D{2} = (LGLWeight1D{1})' .* LGLWeight1D{3};
% xy face
HamDG.grid.LGLWeight2D{3} = (LGLWeight1D{1})' .* LGLWeight1D{2};

% 3D
[X, Y, Z] = ndgrid(LGLWeight1D{1}, LGLWeight1D{2}, LGLWeight1D{3});
HamDG.grid.LGLWeight3D = X .* Y .* Z;


end