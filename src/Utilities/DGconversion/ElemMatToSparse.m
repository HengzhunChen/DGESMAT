function sparseH = ElemMatToSparse(HMat, numElem, sizeHMat)
% ELEMMATTOSPARSE transform an element-wise matrix into a sparse matrix
%
%    sparseH = ElemMatToSparse(HMat, numElem, sizeHMat) convert a
%    element-wise matrix HMat, which is stores in a cell structure, to a
%    sparse matrix whose size is sizeHMat-by-sizeHMat.
%
%    See also HamiltonianDG/ElemMatToSparse.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

numElemTotal = prod(numElem);

numBasisElem = zeros(numElemTotal, 1);
for i = 1 : numElemTotal
    numBasisElem(i) = size(HMat{i, i}, 1);
end

idxI = cell(numElemTotal, numElemTotal);
idxJ = cell(numElemTotal, numElemTotal);
val  = cell(numElemTotal, numElemTotal);

for i = 1 : numElemTotal
    for j = 1 : numElemTotal
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
            idxCols = idxCols(:)';
            value = HMat{i, j}(:)';

            idxI{i, j} = idxRows;
            idxJ{i, j} = idxCols;
            val{i, j} = value;
        end
    end
end

idxI = cell2mat(idxI);
idxJ = cell2mat(idxJ);
val  = cell2mat(val);

sparseH = sparse(idxI, idxJ, val, sizeHMat, sizeHMat);

end