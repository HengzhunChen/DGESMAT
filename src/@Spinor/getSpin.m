function X1 = getSpin(X, component)
% GETSPIN get the different spin components of wavefun
%
%    X1 = getSpin(X, component) returns a Spinor object X1 containing one 
%    component of the Spinor X according to component, which is a number 
%    and should not exceed the number of component of X.
%
%    See also Spinor, Spinor/setSpin.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

componentIdx = [0, cumsum(X.numStateList)];

if component > X.numComponent
    error('wrong input for component of Spinor');
else
    idxSta = componentIdx(component) + 1;
    idxEnd = componentIdx(component + 1);
    X1 = X;
    X1.wavefun = X.wavefun(:, idxSta : idxEnd);
end

X1.numStateList = X.numStateList(component);
X1.numComponent = length(component);

end