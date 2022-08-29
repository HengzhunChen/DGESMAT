function [U, S, V] = svd(X, varargin)
% SPINOR/SVD Singular value decomposition.
%
%    [U, S, V] = svd(X) produces a diagonal matrix S, of the same 
%    dimension as X and with nonnegative diagonal elements in
%    decreasing order, and unitary wave function U and unitary matrix V
%    so that X = U * S * V'.
%  
%    S = svd(X) returns a vector containing the singular values.
%  
%    [U, S, V] = svd(X, 0) produces the "economy size"
%    decomposition. If X is m-by-n with m > n, then only the
%    first n columns of U are computed and S is n-by-n.
%    For m <= n, svd(X, 0) is equivalent to svd(X).
%  
%    [U, S, V] = svd(X, 'econ') also produces the "economy size"
%    decomposition. If X is m-by-n with m >= n, then it is
%    equivalent to svd(X, 0). For m < n, only the first m columns 
%    of V are computed and S is m-by-m.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

if X.trans
    Xpsi = X.wavefun';
else
    Xpsi = X.wavefun;
end

[Umat, S, Vmat] = svd(Xpsi, varargin{:});

if X.trans
    U = Umat;
    V = X;
    V.wavefun = Vmat;
else
    U = X;
    U.wavefun = Umat;
    V = Vmat;
end

end