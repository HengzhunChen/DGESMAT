function HamDG = UpdateDGMatrix(HamDG, vtotLGLDiff)
% HAMILTONIANDG/UPDATEDGMATRIX update the DG Hamiltonian matrix
%
%    HamDG = UpdateDGMatrix(HamDG, vtotLGLDiff) update the DG Hamiltonian 
%    matrix with the same set of adaptive local basis functions, but
%    differenct local potential. This function is used in the innner SCF
%    loop when only the local potential is updated. vtotLGLDiff is the
%    difference of vtot defined on each LGL grid. The contribution of this
%    difference should be added to the DG Hamiltonian matrix.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = HamDG.numElem;

% integration weights
LGLWeight3D = HamDG.grid.LGLWeight3D;

% NOTE: since this is an update process, DO NOT clear the DG Matrix


% *********************************************************************
% Initial setup
% *********************************************************************

% compute global index set
elemBasisIdx = cell(numElem);
elemBasisInvIdx = cell(prod(numElem), 1);
count = 0;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            [~, n] = size(HamDG.basisLGL{i, j, k});
            idxStart = count + 1;
            idxEnd = count + n;
            elemBasisIdx{i, j, k} = idxStart : idxEnd;
            for m = idxStart : idxEnd
                elemBasisInvIdx{m} = [i, j, k];
            end
            count = count + n;
        end
    end
end

HamDG.elemBasisIdx = elemBasisIdx;
HamDG.elemBasisInvIdx = elemBasisInvIdx;
HamDG.sizeHMat = count;



% *********************************************************************
% Update the local potential part
% *********************************************************************

HMat = HamDG.HMat;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            basis = HamDG.basisLGL{i, j, k};
            [~, numBasis] = size(basis);
            
            % skip the calculation is there is no adaptive local basis
            % function
            if numBasis == 0
                continue;
            end
            
            key = i + (j-1) * numElem(1) + (k-1) * numElem(1) * numElem(2);
            localMat = HMat{key, key};
            
            % In all matrix assembly process, only compute the lower
            % triangular matrix and use symmetry later
            
            %
            % local potential part
            %
            vtot = vtotLGLDiff{i, j, k};
            
            tempMat = zeros(numBasis, numBasis);
            for a = 1 : numBasis
                for b = 1 : a
                    tempMat(a, b) = ...
                        sum(basis(:, a) .* LGLWeight3D(:) .* vtot(:) .* basis(:, b));
                end
            end
            
            %
            % symmetrize
            %
            diagVec = diag(tempMat);
            tempMat = tempMat + tempMat';
            tempMat = tempMat - diag(diagVec);

            localMat = localMat + tempMat;
            HMat{key, key} = localMat;
        end
    end
end

HamDG.HMat = HMat;

% HamDG.hasConverted = false;

end