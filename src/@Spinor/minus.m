function Z = minus(X, Y)
% SPINOR/MINUS Operator overload. Minus function for wave function.
%
%    Z = minus(X, Y) returns a Spinor oject with wavefun of X - Y.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

if isa(X, 'Spinor') && isa(Y, 'Spinor')
    Z = X;
    Z.wavefun = X.wavefun - Y.wavefun;
    return;
end

if isa(X, 'Spinor')
    Z = X;
    Z.wavefun = X.wavefun - Y;
    return;
end

if isa(Y, 'Spinor')
    Z = Y;
    Z.wavefun = X - Y.wavefun;
    return;
end

end