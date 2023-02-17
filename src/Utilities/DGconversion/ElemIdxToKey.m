function [i, j, k] = ElemIdxToKey(elemIdx, numElem)
% ELEMIDXTOKEY Index conversion of element.
%
%    [i, j, k] = ElemIdxToKey(elemIdx, numElem) converts the 1-dim index of
%    element elemIdx to 3-dim index [i, j, k] according to number of
%    element in each dimension numElem. Relation of the two types of index
%    are as follow:
%
%    elemIdx = i + (j-1) * numElem(1) + (k-1) * numElem(1)*numElem(2)
%
%    See also ElemKeyToIdx.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


cnt = 0;
for kk = 1 : numElem(3)
    for jj = 1 : numElem(2)
        for ii = 1 : numElem(1)
            cnt = cnt + 1;
            if cnt == elemIdx
                k = kk;
                j = jj;
                i = ii;
                return;
            end
        end
    end
end

end