function [x, w, P, D] = GenerateLGL(N1)
% GENERATELGL generating the unscaled LGL grid points, integration weights, 
%    polynomials and differentiation matrix.
%
%    [x, w, P, D] = GenerateLGL(N1) returns the LGL grid points x, LGL
%    integration weights w, LGL polynomial P, LGL differentiation matrix D 
%    according to LGL grid size N1.
%    NOTE:
%    size(x) = size(w) = size(P, 1|2) = size(D, 1|2) = N1 - 1
%    x, w are saved as column vectors.
%
%    See also LGLMesh, HamiltonianDG, SCFDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

tol = 1e-15;
err = 0;

N = N1 - 1;

x = cos( pi .* (N1-1 : -1 : 0) ./ N );
x = x';  % convert x from row vector to column vector
P = zeros(N1, N1);

iter = 1;
while err > tol || iter == 1
    xold = x;
    
    P(:, 1) = 1;
    P(:, 2) = x;
    for j = 3 : N1
        P(:, j) = ( (2*j-3) .* x .* P(:, j-1) - (j-2) .* P(:, j-2) ) ./ (j-1);
    end
        
    x = xold - ( x .* P(:, N1) - P(:, N1-1) ) ./ (N1 * P(:, N1));
    
    err = max(abs(xold - x));
    
    iter = iter + 1;
end

w = 2 ./ (N * N1 .* P(:, N1).^2 );
w = w';

D = P(:, N1) ./ P(:, N1)' ./ (x - x');
D(1 : (N1+1) : end) = 0;
D(1, 1) = - N * N1 / 4;
D(N1, N1) = N * N1 / 4;

end