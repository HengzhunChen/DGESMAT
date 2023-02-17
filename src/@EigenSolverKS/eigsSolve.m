function EigSol = eigsSolve(EigSol, numEig, eigMaxIter, eigTolerance)
% EIGENSOLVERKS/EIGSSOLVE Solve the eigen-decomposition of planewave-basis 
%    Kohn Sham DFT by function eigs in MATLAB
%
%    EigSol = eigsSolve(EigSol, numEig, eigMaxIter, eigTolerance) computes
%    the partial eigen-decomposition of planewave-basis for Kohn Sham DFT
%    according to number of eigen-pairs to find numEig, maximum number of 
%    iterations eigMaxIter and residual tolerance eigTolerance, and save 
%    the eigen-pairs back to EigSol.eigVal, EigSol.psi.
%
%    See also EigenSolverKS, eigs.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

ntot = EigSol.psi.NumGridTotal();

[V, D, flag] = eigs(@(x) eigsMult(EigSol, x), ...
                    ntot, numEig, 'smallestreal', ...
                    'FailureTreatment', 'keep', ...
                    'IsFunctionSymmetric', true, ...
                    'Tolerance', eigTolerance, ...
                    'MaxIterations', eigMaxIter);

d = real(diag(D));
[sd, id] = sort(d);
V = V(:, id);

if flag ~= 0
    warning('Convergence not reached in eigs!');
end

% save the eigenvalues and eigenvectors back to the EigenSolverKS data
% structure
EigSol.eigVal = sd;
EigSol.psi.wavefun = V;

end


%---------------------------------------------------------------
function y = eigsMult(eigSol, x)

% construct a spinor object from input vector x
X = Spinor(eigSol.psi.domain, 1, x);

% calculate output vector y 
H = eigSol.hamKS;
Y = H * X;
y = Y.wavefun;

end

