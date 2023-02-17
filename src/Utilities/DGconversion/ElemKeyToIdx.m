function elemIdx = ElemKeyToIdx(i, j, k, numElem)
% ELEMKEYTOIDX Index conversion of element.
%
%    elemIdx = ElemKeyToIdx(i, j, k, numElem) converts the 3-dim index of
%    element [i, j, k] to 1-dim index elemIdx according to number of
%    element in each dimension numElem. Relation of the two types of index
%    are as follow:
%
%    elemIdx = i + (j-1) * numElem(1) + (k-1) * numElem(1)*numElem(2)
%
%    See also ElemIdxToKey.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


elemIdx = i + (j-1) * numElem(1) + (k-1) * numElem(1)*numElem(2);

end