function HamDG = Setup(HamDG)
% HAMILTONIANDG/SETUP initializes HamiltonianDG object HamDG with data 
%    from global variable esdfParam.
% 
%    See also HamiltonianDG, Domain, Atom, ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;


% **********************************************************************
% basic parameters
% **********************************************************************

HamDG.domain            = esdfParam.basic.domain;
HamDG.atomList          = esdfParam.basic.atomList;
HamDG.pseudoType        = esdfParam.basic.pseudoType;
HamDG.XCType            = esdfParam.basic.XCType;
HamDG.numExtraState     = esdfParam.basic.numExtraState;
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
HamDG.sizeHMat = prod(HamDG.numElem) * esdfParam.DG.numALBElem(1, 1, 1);
HamDG.hasConvertSparse = 0;



% ***********************************************************************
%  Partition of elements
% ***********************************************************************

HamDG.grid.numLGLGridElem    = esdfParam.DG.numGridLGL;
HamDG.grid.numUniformGridElem = HamDG.domain.numGrid ./ HamDG.numElem;
HamDG.grid.numUniformGridElemFine = HamDG.domain.numGridFine ./ HamDG.numElem;

% Setup the element domains
numElem = HamDG.numElem;
HamDG.domainElem = Domain.empty();
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            key = [i-1, j-1, k-1];
            dm = Domain();
            dm.length = HamDG.domain.length ./ numElem;
            dm.numGrid = HamDG.grid.numUniformGridElem;
            dm.numGridFine = HamDG.grid.numUniformGridElemFine;
            dm.posStart = dm.length .* key;
            HamDG.domainElem(i, j, k) = dm;
        end
    end
end


% Partition by element
HamDG.density     = cell(numElem);
HamDG.densityLGL  = cell(numElem);
HamDG.vext        = cell(numElem);
HamDG.vhart       = cell(numElem);
HamDG.vxc         = cell(numElem);
HamDG.epsxc       = cell(numElem);
HamDG.vtot        = cell(numElem);
HamDG.vtotLGL     = cell(numElem);

% Initialize quantities
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            numGridElemFineTotal = prod(HamDG.grid.numUniformGridElemFine);
            numLGLGridElemTotal  = prod(HamDG.grid.numLGLGridElem);
            
            HamDG.density{i, j, k}     = zeros(numGridElemFineTotal, 1);
            HamDG.densityLGL{i, j, k}  = zeros(numLGLGridElemTotal, 1);
            HamDG.vext{i, j, k}        = zeros(numGridElemFineTotal, 1);
            HamDG.vhart{i, j, k}       = zeros(numGridElemFineTotal, 1);
            HamDG.vxc{i, j, k}         = zeros(numGridElemFineTotal, 1);
            HamDG.epsxc{i, j, k}       = zeros(numGridElemFineTotal, 1);
            HamDG.vtot{i, j, k}        = zeros(numGridElemFineTotal, 1);
            HamDG.vtotLGL{i, j, k}     = zeros(numLGLGridElemTotal, 1);
        end
    end
end

HamDG.gradDensity = cell(dimDef(), 1);
for d = 1 : dimDef()
    HamDG.gradDensity{d} = cell(numElem);
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                HamDG.gradDensity{d}{i, j, k} = zeros(numGridElemFineTotal, 1);
            end
        end
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

HamDG.grid.uniformGridElem = cell(numElem);
HamDG.grid.uniformGridElemFine = cell(numElem);
HamDG.grid.LGLGridElem = cell(numElem);

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            HamDG.grid.uniformGridElem{i, j, k} = ...
                UniformMesh(HamDG.domainElem(i, j, k));
            HamDG.grid.uniformGridElemFine{i, j, k} = ...
                UniformMeshFine(HamDG.domainElem(i, j, k));
            HamDG.grid.LGLGridElem{i, j, k} = ...
                LGLMesh(HamDG.domainElem(i, j, k), HamDG.grid.numLGLGridElem);
        end
    end
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

numLGL = HamDG.grid.numLGLGridElem;
numUniform = HamDG.grid.numUniformGridElem;
numUniformFine = HamDG.grid.numUniformGridElemFine;

EPS = 1e-13;  % small stabilization parameter

for d = 1 : dimDef()
    LGLGrid = HamDG.grid.LGLGridElem{1, 1, 1}{d};
    uniformGrid = HamDG.grid.uniformGridElem{1, 1, 1}{d};
    uniformGridFine = HamDG.grid.uniformGridElemFine{1, 1, 1}{d};
    
    % stabilization constant factor, according to Berrut and Trefethen
    stableFac = 0.25 * HamDG.domainElem(1, 1, 1).length(d);

    % ! Here row vector or column vector is important for computation
        
    lambda = ( LGLGrid - LGLGrid' ) ./ stableFac;
    lambda(1 : numLGL(d)+1 : end) = 1;  % fix diagonal value
    lambda = prod(lambda);
    lambda = 1 ./ lambda;  % size is (1, numLGL(d))
    
    % denominator
    denom = repmat(lambda, numUniform(d), 1) ./ (uniformGrid' - LGLGrid + EPS);
    denom = sum(denom, 2);  % size is (numUniform(d), 1)
    denomFine = repmat(lambda, numUniformFine(d), 1) ./ (uniformGridFine' - LGLGrid + EPS);
    denomFine = sum(denomFine, 2);  % size is (numUniformFine(d), 1)
    
    localMat = (lambda ./ denom) .* (1 ./ (uniformGrid' - LGLGrid + EPS));
    % size is (numUniform(d), numLGL(d))
    localMatFine = (lambda ./ denomFine) .* (1 ./ (uniformGridFine' - LGLGrid + EPS));
    % size is (numUniformFine(d), numLGL(d))
    
    HamDG.grid.LGLToUniformMat{d} = localMat;
    HamDG.grid.LGLToUniformMatFine{d} = localMatFine;
    
% used for reference
%     for i = 1 : numLGL(d)
%         lambda(i) = 1.0;
%         for j = 1 : numLGL(d)
%             if j ~= i 
%                 lambda(i) = lambda(i) * (LGLGrid(i) - LGLGrid(j)) / stableFac; 
%             end
%         end
%         lambda(i) = 1.0 / lambda(i);
%         for j = 1 : numUniform(d)
%             denom(j) = denom(j) + lambda(i) / ( uniformGrid(j) - LGLGrid(i) + EPS );
%         end
%         for j = 1 : numUniformFine(d)
%             denomFine(j) = denomFine(j) + lambda(i) / ( uniformGridFine(j) - LGLGrid(i) + EPS );
%         end
%     end
% 
%     for i = 1 : numLGL(d)
%         for j = 1 : numUniform(d)
%             localMat( j, i ) = (lambda(i) / ( uniformGrid(j) - LGLGrid(i) + EPS )) / denom(j); 
%         end
%         for j = 1 : numUniformFine(d)
%             localMatFine( j, i ) = (lambda(i) / ( uniformGridFine(j) - LGLGrid(i) + EPS )) / denomFine(j);
%         end
%     end
    
end

% ----------------------------------------------------------------------
% Compute the LGL weights at 1D, 2D, 3D
% ----------------------------------------------------------------------

dmElemlength = HamDG.domainElem(1, 1, 1).length;
numGrid = HamDG.grid.numLGLGridElem;

% compute the integration weights
% 1D
for d = 1 : dimDef()
    [~, LGLWeight1D, ~, ~] = GenerateLGL(numGrid(d));
    HamDG.grid.LGLWeight1D{d} = 0.5 * dmElemlength(d) .* LGLWeight1D;
end

% 2D: faces labeled by normal vectors, i.e.
% yz face : 1
% xz face : 2
% xy face : 3

% yz face
HamDG.grid.LGLWeight2D{1} = (HamDG.grid.LGLWeight1D{2})' .* HamDG.grid.LGLWeight1D{3}; 
% xz face
HamDG.grid.LGLWeight2D{2} = (HamDG.grid.LGLWeight1D{1})' .* HamDG.grid.LGLWeight1D{3};
% xy face
HamDG.grid.LGLWeight2D{3} = (HamDG.grid.LGLWeight1D{1})' .* HamDG.grid.LGLWeight1D{2};

% 3D
[X, Y, Z] = ndgrid(HamDG.grid.LGLWeight1D{1}, HamDG.grid.LGLWeight1D{2}, HamDG.grid.LGLWeight1D{3});
HamDG.grid.LGLWeight3D = X .* Y .* Z;


end