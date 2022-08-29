function Y = uminus(X)
% SPINOR/UMINUS Operator overload. Uminus function for wave function.
%
%    Y = uminus(X) returns a Spinor with wave function as -X.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.wavefun = -X.wavefun;

end