function Y = power(X, p)
% SPINOR/POWER Operator overload. Power function for wave function
%
%    Y = power(X, p) returns a Spinor with wave function of X.^p.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.wavefun = (X.wavefun) .^ p;

end