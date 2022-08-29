function globalVec = ElemVecToGlobal(elemVec, numGrid, numElem)
% ELEMVECTOGLOBAL transform a vector in element-wise structure to vector 
%    over global domain.
%
%    globalVec = ElemVecToGlobal(elemVec, numGrid, numElem) convert a 
%    element-wise vector elemVec, which is saved as a cell, to a vector 
%    over global domain globalVec. numGrid could be coarse or fine grid 
%    over global domain.
%
%    See also GlobalVecToElem.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

numGridElem = numGrid ./ numElem;

elemTns = cell(numElem);

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            elemTns{i, j, k} = reshape(elemVec{i, j, k}, numGridElem);
        end
    end
end

globalTns = cell2mat(elemTns);

globalVec = reshape(globalTns, [numGrid(1)*numGrid(2)*numGrid(3), 1]);

end