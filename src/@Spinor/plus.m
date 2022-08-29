function Z = plus(X, Y)
% SPINOR/PLUS Operator overload. Plus function for wave function
%
%    Z = plus(X, Y) returns a Spinor with wave function of X + Y.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

if isa(X, 'Spinor') && isa(Y, 'Spinor')
    Z = X;
    Z.wavefun = X.wavefun + Y.wavefun;
    return;
end

if isa(X, 'Spinor')
    Z = X;
    Z.wavefun = X.wavefun + Y;
    return;
end

if isa(Y, 'Spinor')
    Z = Y;
    Z.wavefun = X + Y.wavefun;
    return;
end

end