function elemVec = GlobalVecToElem(globalVec, numGrid, numElem)
% GLOBALVECTOELEM transform vector over global domain to element-wise
%
%    elemVec = GlobalVecToElem(globalVec, numGrid, numElem) convert a 
%    vector over global domain globalVec into the element-wise vector 
%    structure, which is stored as a cell elemVec. numGrid could be coarse 
%    or fine grid over global domain.
%
%    See also ElemVecToGlobal.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

numGridElem = numGrid ./ numElem;

globalTns = reshape(globalVec, numGrid);

dim1Dist = numGridElem(1) * ones(numElem(1), 1);
dim2Dist = numGridElem(2) * ones(numElem(2), 1);
dim3Dist = numGridElem(3) * ones(numElem(3), 1);
eleTns = mat2cell(globalTns, dim1Dist, dim2Dist, dim3Dist);

elemVec = cell(numElem);

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            elemVec{i, j, k} = reshape(eleTns{i, j, k}, ...
                [numGridElem(3)*numGridElem(2)*numGridElem(1), 1]);
        end
    end
end

end