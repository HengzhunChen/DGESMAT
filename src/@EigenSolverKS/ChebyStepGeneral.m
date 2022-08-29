function EigSol = ChebyStepGeneral(EigSol, numEig, filterOrder)
% EIGENSOLVERKS/CHEBYSTEPGENERAL General Chebyshev filter step of Chebyshev 
%    filtering subspace iteration (CheFSI) for eigen-decomposition of 
%    planewave-basis Kohn Sham DFT
%
%    EigSol = ChebyStepGeneral(EigSol, numEig, filterOrder) performs the 
%    general step of Chebyshev subspace iteration (CheFSI) of 
%    planewave-basis for Kohn Sham DFT according to number of eigen-pairs 
%    to find numEig, filter order filterOrder, and save the eigen-pairs 
%    and residual value back to EigSol.eigVal, EigSol.psi, EigSol.resVal.
%
%    See also EigenSolverKS, EigenSolverKS/ChebyStepFirst.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


InfoPrint(0, "\n Subsequent CheFSI for PWDFT ... \n");

% ***********************************************************************
% Initialization
% ***********************************************************************

width = EigSol.psi.NumStateTotal();

if numEig > width
    msg = "Number of eigenvalues requested = " + num2str(numEig) + ...
          ", which is larger than the number of columns in psi = " + ...
          num2str(width);
    error(msg);
end

% TODO: in what follows, we do not consider multiple spinor components
% For now, just a spin unpolarized setup is considered

% TODO: in what follows, we assume real value wavefunctions

%
% Step 1: Obtain the upper bound and the Ritz values (for the lower bound)
% using the Lanczos estimator
%

InfoPrint(0, '\n Estimating the upper bound ...');
numLanczosSteps = 5;
[b_up, ~] = ChebyBound(EigSol, numLanczosSteps);


%
% Step 2: Set up the lower bound from previous Rayleigh-Ritz / Eigensolver
% values
%
a_L = EigSol.eigVal(1);  % scale factor to prevent potential overflow 
b_low = EigSol.eigVal(width);

InfoPrint(0, "\n Upper bound (to be mapped to +1) = %1.8e", b_up);
InfoPrint(0, "\n Lower bound (to be mapped to -1) = %1.8e", b_low);
InfoPrint(0, "\n Lowest eigenvalue = %1.8e", a_L);


% ***********************************************************************
% Step 3: Main loop
% ***********************************************************************

X = EigSol.psi;

% 
% Step 3a: Compute the filtered block of vectors
% This always works on the vector of psi
%
% option 1:
X = ChebyFilterScaled(EigSol, X, filterOrder, b_low, b_up, a_L);
% option 2:
% X = ChebyFilter(EigSol, X, filterOrder, b_low, b_up);


%
% Step 3b: orthonormalize
%
% Option 1: Cholesky factorization 
% XTX = X' * X;
% XTX = (XTX + XTX') / 2;
% U = chol(XTX); 
% X = X / U;
% option 2: QR decomposition
[X, ~] = qr(X, 0);


%
% Step 3c: Rayleigh-Ritz step
%
H = EigSol.hamKS;

HX = H * X;
XTHX = X' * HX;
XTHX = (XTHX + XTHX') / 2;
[V, D] = eigs(XTHX, width);
d = real(diag(D));
[eigVals, id] = sort(d);
Q = V(:, id);


%
% step 3d: Subspace rotation step: psi <-- psi * Q
% where Q are eigenvectors from the Rayleigh-Ritz step
%
X = X * Q;


% **********************************************************************
% Post processing
% **********************************************************************

% save data back to EigenSolverKS structure
EigSol.eigVal = eigVals;
EigSol.psi = X;

% compute the residual of the eigenvectors
Res = HX - X * (X' * H * X); 

% compute the norm of residuals
resNorm = vecnorm(Res, 2, 1);     
resMax = max(resNorm);
resMin = min(resNorm);

InfoPrint(0, '\n Maximum residual = %1.8e', resMax);
InfoPrint(0, '\n Minimum residual = %1.8e \n\n', resMin);

EigSol.resVal = resNorm;


end