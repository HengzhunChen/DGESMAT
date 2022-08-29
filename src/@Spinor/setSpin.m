function X = setSpin(X, component, X1)
% SETSPIN set the different spin components of wavefun
%
%    X = setSpin(X, component, X1) delivers the wavefun in Spinor X1 to one
%    component of wavefun in Spinor X, where component is a number and
%    should not exceed the number of component of X.
%
%    See also Spinor, Spinor/getSpin.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

componentIdx = [0, cumsum(X.numStateList)];

if component > X.numComponent
    error('wrong input for component of Spinor');
else
    idxSta = componentIdx(component) + 1;
    idxEnd = componentIdx(component + 1);
    X.wavefun = ...
        [X.wavefun(:, 1 : idxSta-1), X1.wavefun, X.wavefun(:, idxEnd+1 : end)];
end

X.numStateList(component) = sum(X1.numStateList);

end