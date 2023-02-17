function globalVec = ElemVecToGlobal(elemVec, numGrid, numElem)
% ELEMVECTOGLOBAL transform a vector in element-wise structure to vector 
%    over global domain.
%
%    globalVec = ElemVecToGlobal(elemVec, numGrid, numElem) convert a 
%    element-wise vector elemVec, which is saved in a 1-dim cell with 
%    length numElemTotal, to a vector over global domain globalVec. 
%    numGrid could be coarse or fine grid over global domain.
%
%    NOTE: for input elemVec, it can also have cell member with a matrix 
%    whose columns are all vectors over corresponding element. 
%
%    See also GlobalVecToElem.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numGridTotal = prod(numGrid);
numCol = size(elemVec{1}, 2);

globalVec = zeros(numGridTotal, numCol);
elemVec = reshape(elemVec, numElem);

numGridElem = numGrid ./ numElem;

for n = 1 : numCol
    elemTns = cell(numElem);
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                elemTns{i, j, k} = reshape(elemVec{i, j, k}(:, n), numGridElem);
            end
        end
    end
    globalTns = cell2mat(elemTns);
    globalVec(:, n) = reshape(globalTns, [numGrid(1)*numGrid(2)*numGrid(3), 1]);
end

end