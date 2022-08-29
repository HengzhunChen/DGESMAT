function Y = ChebyFilter(EigSol, X, m, a, b)
% EIGENSOLVERKS/CHEBYFILTER subroutine to perform the unscaled Chebyshev
%    filtering
%
%    Y = ChebyFilter(EigSol, X, m, a, b) perform the unscaled Chebyshev
%    filter to X according to filter order m and eigenvalues interval 
%    [a, b].
%
%    See also EigenSolverKS, EigenSolveKS/ChebyFilterScaled.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

F = EigSol.fft;
H = EigSol.hamKS;

e = (b - a) / 2;
c = (a + b) / 2;

% apply energy cutoff filter since we are starting from a random guess
if H.ecutFilter.isApplyFilter == 1
    Y = F * X;
    idx = (F.gkk ./ 2) > H.ecutFilter.wfnCutoff;
    Y(idx, :) = 0;
    X = F' * Y;
    X = real(X);
end

Y = (H * X - c * X) / e;

for i = 2 : m
    Yt = (H * Y - c * Y) * (2 / e) - X;
    X = Y;
    Y = Yt;
end

end