function Z = mtimes(X, Y)
% SPINOR/MTIMES Operator overload. Mtimes function for wave function.
%
%    M = mtimes(X, Y) returns a matrix as the multiplication of two wave 
%    functions if Y is a Spinor object.
%
%    Z = mtimes(M, Y) or Z = mtimes(X, M) returns a wave function as the
%    multiplication of a matrix or scalar with a wave function.
%
%    Z = mtimes(Xtrans, H) returns a Spinor X' * H = (H * X)', where H is
%    a HamiltonianKS object and Xtran = X'.
%
%    NOTE: the left multiplication of a non-scalar matrix with a wave
%    function returns a regular matrix.
%
%    See also Spinor, HamiltonianKS.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


if isa(X, 'Spinor') && isa(Y, 'Spinor')
    if X.trans
        Xpsi = X.wavefun';
    else
        Xpsi = X.wavefun;
    end
    if Y.trans
        Ypsi = Y.wavefun';
    else
        Ypsi = Y.wavefun;
    end
    Z = Xpsi * Ypsi;
    return;
end

if isa(X, 'Spinor') && isa(Y, 'HamiltonianKS')
    if X.trans
        Z = Y * X';
        Z = Z';
    else
        error("syntax X * H error, it should be X' * H"); 
        % here X is a Spinor object and H is a HamiltonianKS object
    end
    return;
end

if isa(X, 'Spinor')
    if X.trans
        Z = X.wavefun' * Y;
    else
        Z = X;
        Z.wavefun = X.wavefun * Y;
    end
    return;
end

if isa(Y, 'Spinor')
    if Y.trans
       Z = Y;
       Z.wavefun = Y.wavefun * X';
    else
        if isscalar(X)
            Z = Y;
            Z.wavefun = X * Y.wavefun;
        else
            Z = X * Y.wavefun;
        end
    end
    return;
end

end