function [b_up, ritzValues] = ChebyBound(EigSol, numLanczosSteps)
% EIGENSOLVERKS/CHEBYBOUND subroutine to compute the bound of eigenvalues
%    used in Chebyshev filter.
%
%    [b_up, ritzValues] = ChebyBound(EigSol, numLanczosSteps) computes the
%    upper bound of eigenvalues b_up and ritz values ritzValues which is
%    used to esimate the lower bound of eigenvalues.
%
%    See also EigenSolverKS, EigenSolverKS/ChebyStepFirst,
%    EigenSolverKS/ChebyStepGeneral.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

F = EigSol.fft;
H = EigSol.hamKS;

% Setup a temporary Spinor with random numbers
numState = 1;
v = Spinor(F.domain, numState, 0.0);
v.wavefun = rand( size(v.wavefun) );

% Normalize the temporary spinor
v = v ./ norm(v, 2);

% apply energy cutoff filter since we are starting from random guess
if H.ecutFilter.isApplyFilter == 1
    y = F * v;
    idx = (F.gkk ./ 2) > H.ecutFilter.wfnCutoff;
    y(idx, :) = 0;
    v = F' * y;
    v = real(v);
end

f = H * v;
alpha = f' * v;
f = f - alpha * v;

T = zeros(numLanczosSteps, numLanczosSteps);
T(1, 1) = alpha;

for j = 2 : numLanczosSteps
    beta = norm(f, 2);
    v0 = v;
    v = f ./ beta;
    f = H * v;
    f = f - beta * v0;
    alpha = f' * v;
    f = f - alpha * v;
    T(j, j-1) = beta;
    T(j-1, j) = beta;
    T(j, j) = alpha;
end

% Solve the eigenvalue problem for the Ritz values
T = (T + T') / 2;
[~, eigval] = eigs(T, numLanczosSteps);
d = real(diag(eigval));
[ritzValues, ~] = sort(d);

b_up = ritzValues(numLanczosSteps) + norm(f, 2);
% b_up = norm(T, 2) + norm(f, 2);


end