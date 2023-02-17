function X = mldivide(H, B)
% HAMILTONIANKS/MLDIVIDE Operator overload solving HX = B, where H is a
%    HamiltonianKS object and both X and B are Spinor objects.
%
%    Usage: X = H \ B;
%    B is a Spinor object containing wave function info.
%
%    See also HamiltonianKS, Spinor, HamiltonianKS/mtimes.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

b = B.wavefun;
x = zeros(size(b));

for i = 1 : ncols(B)
    [x(:, i), flag] = gmres( @(p) mult(H, p), b(:, i), [], 1e-6, 300 );

    if flag ~= 0
        warning('Convergence not reached in gmres');
    end
end

X = B;
X.wavefun = x;

end


% Auxilary function
function y = mult(H, x)

X = Spinor(H.domain, 1, 1, x);
Y = H * X;
y = Y.wavefun;

end