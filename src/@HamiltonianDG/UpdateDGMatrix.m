function HamDG = UpdateDGMatrix(HamDG, vtotLGLDiff)
% HAMILTONIANDG/UPDATEDGMATRIX update the DG Hamiltonian matrix
%
%    HamDG = UpdateDGMatrix(HamDG, vtotLGLDiff) update the DG Hamiltonian 
%    matrix with the same set of adaptive local basis functions, but
%    differenct local potential. This function is used in the inner SCF
%    loop when only the local potential is updated. vtotLGLDiff is the
%    difference of vtot defined on each LGL grid. The contribution of this
%    difference should be added to the DG Hamiltonian matrix.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = HamDG.numElem;
numElemTotal = prod(numElem);

% integration weights
LGLWeight3D = HamDG.grid.LGLWeight3D;
LGLWeight3D = reshape(LGLWeight3D, [], 1);

% NOTE: since this is an update process, DO NOT clear the DG Matrix


% *********************************************************************
% Initial setup
% *********************************************************************

% compute global index set
elemBasisIdx = cell(numElemTotal, 1);
elemBasisInvIdx = cell(numElemTotal, 1);
count = 0;

for elemIdx = 1 : numElemTotal
    [~, n] = size(HamDG.basisLGL{elemIdx});
    idxStart = count + 1;
    idxEnd = count + n;
    elemBasisIdx{elemIdx} = idxStart : idxEnd;
    for m = idxStart : idxEnd
        elemBasisInvIdx{m} = elemIdx;
    end
    count = count + n;
end

HamDG.elemBasisIdx = elemBasisIdx;
HamDG.elemBasisInvIdx = elemBasisInvIdx;
HamDG.sizeHMat = count;


% *********************************************************************
% Update the local potential part
% *********************************************************************

HMat = HamDG.HMat;
basisLGL = HamDG.basisLGL;

for elemIdx = 1 : numElemTotal
    basis = basisLGL{elemIdx};
    [~, numBasis] = size(basis);
    
    % skip the calculation is there is no adaptive local basis
    % function
    if numBasis == 0
        continue;
    end
    
    localMat = HMat{elemIdx, elemIdx};
        
    %
    % local potential part
    %
    vtot = vtotLGLDiff{elemIdx};
    localMat = localMat + basis' * ((LGLWeight3D .* vtot) .* basis); 
    
    HMat{elemIdx, elemIdx} = localMat;
end

HamDG.HMat = HMat;


end