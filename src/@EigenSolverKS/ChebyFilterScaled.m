function Y = ChebyFilterScaled(EigSol, X, m, a, b, a_L)
% EIGENSOLVERKS/CHEBYFILTERSCALED subroutine to perform the scaled 
%    Chebyshev filtering
%
%    Y = ChebyFilterScaled(EigSol, X, m, a, b, a_L) perform the scaled 
%    Chebyshev filter to X according to filter order m and eigenvalues 
%    interval [a, b], a_L is scaling factor smaller than interval endpoint 
%    a to prevent potential overflow.
%
%    See also EigenSolverKS, EigenSolveKS/ChebyFilter.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

F = EigSol.fft;
H = EigSol.hamKS;

e = (b - a) / 2;
c = (a + b) / 2;
sigma = e / (c- a_L);
tau = 2.0 / sigma;

% apply energy cutoff filter since we are starting from a random guess
if H.ecutFilter.isApplyFilter == 1
    Y = F * X;
    idx = (F.gkk ./ 2) > H.ecutFilter.wfnCutoff;
    Y(idx, :) = 0;
    X = F' * Y;
    X = real(X);
end

Y = (H * X - c * X) .* (sigma / e);

for i = 2 : m
    sigma_new = 1 / (tau - sigma);
    Yt = (H * Y - c * Y) .* (2 * sigma_new / e) - (sigma * sigma_new) * X;
    X = Y;
    Y = Yt;
    sigma = sigma_new;
end

end