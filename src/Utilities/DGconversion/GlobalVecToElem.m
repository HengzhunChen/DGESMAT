function elemVec = GlobalVecToElem(globalVec, numGrid, numElem)
% GLOBALVECTOELEM transform vector over global domain to element-wise
%
%    elemVec = GlobalVecToElem(globalVec, numGrid, numElem) convert a 
%    vector over global domain globalVec into the element-wise vector 
%    structure, which is stored as a 1-dim cell elemVec with length 
%    numElemTotal. numGrid could be coarse or fine grid over global 
%    domain.
%
%    NOTE: input globalVec can also be a matrix with each column a 
%    vector over global domain. 
%
%    See also ElemVecToGlobal.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

numGridElem = numGrid ./ numElem;

dim1Dist = numGridElem(1) * ones(numElem(1), 1);
dim2Dist = numGridElem(2) * ones(numElem(2), 1);
dim3Dist = numGridElem(3) * ones(numElem(3), 1);

elemVec = cell(numElem);

for n = 1 : size(globalVec, 2)
    globalTns = reshape(globalVec(:, n), numGrid);
    eleTns = mat2cell(globalTns, dim1Dist, dim2Dist, dim3Dist);
    
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                elemVec{i, j, k}(:, n) = reshape(eleTns{i, j, k}, ...
                    [numGridElem(3)*numGridElem(2)*numGridElem(1), 1]);
            end
        end
    end
end

elemVec = reshape(elemVec, [], 1);

end