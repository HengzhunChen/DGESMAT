function EigSol = ChebyStepFirst(EigSol, numEig, eigMaxiter, filterOrder)
% EIGENSOLVERKS/CHEBYSTEPFIRST First Chebyshev filter step of Chebyshev 
%    filtering subspace iteration (CheFSI) for eigen-decomposition of 
%    planewave-basis Kohn Sham DFT
%
%    EigSol = ChebyStepFirst(EigSol, numEig, eigMaxiter, filterOrder) 
%    performs the first step of Chebyshev subspace iteration (CheFSI) of 
%    planewave-basis for Kohn Sham DFT according to number of eigen-pairs 
%    to find numEig, maximum number of iteration eigMaxIter, filter order
%    filterOrder, and save the eigen-pairs and residual value back to 
%    EigSol.eigVal, EigSol.psi, EigSol.resVal.
%
%    See also EigenSolverKS, EigenSolverKS/ChebyStepGeneral.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


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

InfoPrint(0, '\n Estimating the upper bound ... \n');
numLanczosSteps = 6;
[b_up, RitzValues] = ChebyBound(EigSol, numLanczosSteps);


%
% Step 2: Set up the lower bound and the filter scale
%
a_L = RitzValues(1);  % a_L is the minimal one of Ritz values 
beta = 0.5;  % 0.5 <= beta < 1.0
b_low = beta * RitzValues(1) + (1.0 - beta) * RitzValues(numLanczosSteps);


% ***********************************************************************
% Step 3: Main loop
% ***********************************************************************

iterMax = eigMaxiter;

for iter = 1 : iterMax
    InfoPrint(0, ' First CheFSI for PWDFT cycle %d of %d \n', iter, iterMax);
    
    InfoPrint(0, " Upper bound (to be mapped to +1) = ", b_up);
    InfoPrint(0, " Lower bound (to be mapped to -1) = ", b_low);
    InfoPrint(0, " Lowest eigenvalue = ", a_L);

    X = EigSol.psi;

    % 
    % Step 3a: Compute the filtered block of vectors
    % This always works on the vectors in psi
    %
    X = ChebyFilterScaled(EigSol, X, filterOrder, b_low, b_up, a_L);

    % 
    % Step 3b: Orthonormalize
    %
    % Option 1: Cholesky factorization 
%     XTX = X' * X;
%     XTX = (XTX + XTX') / 2;
%     U = chol(XTX); 
%     X = X / U;
    % option 2: QR factorization
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
    % Step 3d: Subspace rotation step: psi <-- psi * Q
    % where Q are the eigenvectors from the Rayleigh-Ritz step
    %
    X = X * Q;
    EigSol.psi = X;

    %
    % Step 3e: Reset the upper and lower bounds using the results of the
    % Rayleigh-Ritz step
    %
    b_low = eigVals(width);
    a_L = eigVals(1);


    % *******************************************************************
    % Post processing
    % *******************************************************************

    % compute the residual of the eigenvectors
    Res = HX - X * (X' * H * X);
    resNorm = vecnorm(Res, 2, 1);
    resMax = max(resNorm);
    resMin = min(resNorm);
    
    InfoPrint(0, ' Maximum residual = %1.8e ', resMax);
    InfoPrint(0, '\n Minimum residual = %1.8e \n\n', resMin);
    
    EigSol.resVal = resNorm;
end

% save the eigenvalues to the eigensolver data structure
EigSol.eigVal = eigVals;

end