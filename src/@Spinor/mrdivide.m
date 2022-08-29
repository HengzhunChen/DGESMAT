function Z = mrdivide(X, R)
% SPINOR/MRDIVIDE Operator overload. Mrdivide function for wave function. 
%
%    Z = mrdivide(X, R) returns a Spinor object with wave function of X/R.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Z = X;
Z.wavefun = X.wavefun / R;

end