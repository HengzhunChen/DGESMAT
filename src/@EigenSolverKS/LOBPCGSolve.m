function EigSol = LOBPCGSolve(EigSol, numEig, eigMaxIter, eigMinTolerance, eigTolerance)
% EIGENSOLVERKS/LOBPCGSOLVE LOBPCG solver of eigen-decomposition for 
%    planewave-basis Kohn Sham DFT
%
%    EigSol = LOBPCGSolve(EigSol, numEig, eigMaxIter, eigMinTolerance, 
%    eigTolerance) computes the partial eigen-decomposition of 
%    planewave-basis for Kohn Sham DFT according to number of eigen-pairs 
%    to find numEig, maximum number of iteration eigMaxIter, minimum 
%    tolerance must be reached eigMinTolerance and residual tolerance 
%    eigTolerance, and save the eigen-pairs and residual value back to 
%    EigSol.eigVal, EigSol.psi, EigSol.resVal.
%
%    See also EigenSolverKS.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


timeStart = tic;

% ************************************************************************
% Initialization
% ************************************************************************

width = EigSol.psi.NumStateTotal();

if numEig > width
    msg = "Number of eigenvalues requested = " + num2str(numEig) + ...
          ", which is larger than the number of columns in psi = " + num2str(width);
    error(msg);
end

X = EigSol.psi;
H = EigSol.hamKS;
F = EigSol.fft;


% S = ( X | W | P ) is a triplet used for LOBPCG.
% W is the preconditioned residual
% intermediate matrix temporarily used
% AMat = S' * H * S = S' * (AS)
% BMat = S' * S

isRestart = false;
isConverged = false;

% numLocked is the number of converged vectors
% NOTE: never perform locking in this version
numLocked = 0;


% ************************************************************************
% Main loop
% ************************************************************************

% ----------------------- Orthogonalization ---------------------------

% Option 1: Cholesky factorization 
XTX = X' * X;
U = chol(XTX); 
X = X / U;

% Option 2: QR decomposition
% [X, ~] = qr(X, 0);

AX = H * X;


% --------------------- Start the main loop ----------------------------
iter = 0;
InfoPrint(0, '\nMinmum tolerance is %1.8e \n', eigMinTolerance);

P = zeros(nrows(X), ncols(X));

resMin = eigMinTolerance + 10;

while iter < (10 * eigMaxIter) && ( iter < eigMaxIter || resMin > eigMinTolerance )
    iter = iter + 1;
    
    if iter == 1 || isRestart
        numSet = 2;  % Steepest descent (Davidson), only use (X | W)
    else
        numSet = 3;  % Conjugate gradient, use all the triplet (X | W | P)
    end
    
    % compute the residual
    XTHX = X' * AX;
    R = AX - X * XTHX;
        
    % compute the norm of the residual
    resNorm = vecnorm(R, 2, 1);
    resNorm = resNorm ./ max(1, abs(diag(XTHX))');     
    resMax = max(resNorm);
    resMin = min(resNorm);
        
    if resMax < eigTolerance
        isConverged = true;
        break;
    end    
    
    % compute the preconditioned residual W = T * R
    W = F' * (F*R .* F.TeterPrecond);
    W = real(W);
    % normalize the preconditioned residual
    normW = vecnorm(W, 2, 1);
    W = W ./ normW;
    
    % normalize the conjugate direction
    if numSet == 3
        normP = vecnorm(P, 2, 1);
        P = P ./ normP;
        AP = AP ./ normP;
    end
    
    % compute AMat = S'* H * S = S' * AS
    AW = H * W;
    
    % separate the save run time
    if numSet == 2
        XTHX = X' * AX;
        XTHW = X' * AW;
        WTHW = W' * AW;
        AMat = [XTHX,  XTHW; ...
                XTHW', WTHW];
    end
    
    if numSet == 3
        XTHX = X' * AX;
        XTHW = X' * AW;
        XTHP = X' * AP;
        WTHP = W' * AP;
        WTHW = W' * AW;
        PTHP = P' * AP;
        AMat = [XTHX,  XTHW,  XTHP; ...
                XTHW', WTHW,  WTHP; ...
                XTHP', WTHP', PTHP];
    end
    
    % Compute BMat (overlap matrix)
    if numSet == 2
        S = [X, W];
        BMat = S' * S;
    end
    if numSet == 3
        S = [X, W, P];
        BMat = S' * S;
    end
    
    
    % Rayleigh-Ritz procedure
    % AMat * C = BMat * C * Lambda
    % Assuming the dimension (needed) for C is width * width, then
    %     ( C_X )
    %     ( --- )
    % C = ( C_W )
    %     ( --- )
    %     ( C_P )
    
    numActive = width - numLocked;    

    if numSet == 3
        % conjugate gradient
        numCol = width + 2 * numActive;
    else
        numCol = width + numActive;
    end
    
        
    % Symmetrize A and B first. This is important.
    AMat = (AMat + AMat') / 2;
    BMat = (BMat + BMat') / 2;
    
    [V, D] = eig(BMat);    
    d = real(diag(D));
    [sigma2, id] = sort(d); % in ascending order
    sigma2 = sigma2(1 : numCol);
    id = id(1 : numCol);
    BMat = V(:, id);

    
    idx = find( sigma2 ./ sigma2(end) > 1e-8 );
    numKeep = length(sigma2) - idx(1) + 1;
    
    if numKeep < width
        error('There are not enough number of columns.');
    end
    
    invsigma = 1.0 ./ sqrt( sigma2(end-numKeep+1 : end) );
    
    
    % evaluate S^{-1/2} (U^T A U) S^{-1/2}
    BMat = BMat(:, end-numKeep+1 : end);
    AMat1 = BMat' * AMat * BMat;  % size (numKeep, numKeep)
    AMat1 = AMat1 .* (invsigma * invsigma');
    
    % Solve the standard eigenvalue problem
    AMat1 = (AMat1 + AMat1') / 2;
    [V, eigval] = eig(AMat1);
    d = real(diag(eigval));
    [~, id] = sort(d);
    id = id(1 : numKeep);
    AMat1 = V(:, id);
    
    % Compute the correct eigenvactors and save them in AMat
    AMat1 = AMat1 .* invsigma;    
    AMat1 = BMat * AMat1;    
    AMat(:, 1:numKeep) = AMat1;
    
    
    if numSet == 2
        CX = AMat(          1  :    width,  1 : width );
        CW = AMat( (  width+1) : (2*width), 1 : width );
    end
    if numSet == 3
        CX = AMat(          1  :    width,  1 : width );
        CW = AMat( (  width+1) : (2*width), 1 : width );
        CP = AMat( (2*width+1) : (3*width), 1 : width );
    end
    
    if numSet == 2
        % update the eigenvectors
        X = X * CX + W * CW;        
        P = W;
    end
    if numSet == 3
        % compute the conjugate direction
        P = W * CW + P * CP;        
        % update the eigenvectors
        X = X * CX + P;
    end
    
    % update AX and AP
    if numSet == 2
        AX = AX * CX + AW * CW;        
        AP = AW;
    end    
    if numSet == 3
        AP = AW * CW + AP * CP;
        AX = AX * CX + AP;
    end            
end


% ************************************************************************
% Post processing
% ************************************************************************

% Obtain the eigenvalues and eigenvectors
% XTHX sholud now contain the matrix X' * (AX), and X is an orthonormal set

XTHX = (XTHX + XTHX') / 2;
[V, D] = eig(XTHX);
d = real(diag(D));
[eigValS, id] = sort(d);
XTHX = V(:, id);

X = X * XTHX;

% Save the eigenvalues and eigenvectors back to the EigenSolverKS data
% structure
EigSol.eigVal = eigValS;
EigSol.resVal = resNorm;
EigSol.psi = X;

timeEnd = toc(timeStart);

if isConverged
    InfoPrint(0, "After %d iterations, LOBPCG has converged. \n", iter);
    InfoPrint(0, "The maximum norm of the residual is %1.8e \n", resMax); 
    InfoPrint(0, "The minimum norm of the residual is %1.8e \n", resMin);
    InfoPrint(0, "Time for this LOBPCG call is %8f [s] \n\n", timeEnd);
else
    InfoPrint(0, "After %d iterations, LOBPCG did not converge. \n", iter);
    InfoPrint(0, "The maximum norm of the residual is %1.8e \n", resMax); 
    InfoPrint(0, "The minimum norm of the residual is %1.8e \n", resMin);
    InfoPrint(0, "Time for this LOBPCG call is %8f [s] \n\n", timeEnd);    
end
   

end