function HX = mtimes(HamDG, X)
% HAMILTONIANDG/MTIMES overload mutiplication operator between 
%    HamiltonianDG object H and matrix X.  
%
%    Usage: HX = H * X;
%    H is a HamiltonianDG and X is a matrix.
%
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% NOTE: temporarily implement, may be fixed later. 

HMat = HamDG.HMat;
numElem = HamDG.numElem;
numElemTotal = prod(numElem);

numBasisElem = zeros(numElemTotal, 1);
for i = 1 : numElemTotal
    numBasisElem(i) = size(HMat{i, i}, 1);
end


sparseH = sparse(numElemTotal, numElemTotal);

for i = 1 : prod(numElem)
    for j = 1 : prod(numElem)
        if ~isempty(HMat{i, j})
            numRows = numBasisElem(i);
            numCols = numBasisElem(j);
            
            idxRowStart = sum(numBasisElem(1 : (i-1))) + 1;
            idxRowEnd = idxRowStart + numRows - 1;
            idxRow = idxRowStart : idxRowEnd;
            
            idxColStart = sum(numBasisElem(1 : (j-1))) + 1;
            idxColEnd = idxColStart + numCols - 1;
            idxCol = idxColStart : idxColEnd;
            
            idxRows = repmat(idxRow, 1, numCols);
            idxCols = repmat(idxCol, numRows, 1);
            idxCols = idxCols(:);
            value = HMat{i, j}(:);
            tempM = sparse(idxRows, idxCols, value, numElemTotal, numElemTotal);
            sparseH = sparseH + tempM;
        end
    end
end

HX = sparseH * X;


end