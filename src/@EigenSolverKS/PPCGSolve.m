function EigSol = PPCGSolve(EigSol, numEig, eigMaxIter)
% EIGENSOLVERKS/PPCGSOLVE PPCG solver of eigen-decomposition for 
%    planewave-basis Kohn Sham DFT
%
%    EigSol = PPCGSolve(EigSol, numEig, eigMaxIter) computes the partial 
%    eigen-decomposition of planewave-basis for Kohn Sham DFT according to 
%    number of eigen-pairs to find numEig, maximum number of iteration 
%    eigMaxIter, and save the eigen-pairs and residual value back to 
%    EigSol.eigVal, EigSol.psi, EigSol.resVal.
%
%    See also EigenSolverKS.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

sbSize = esdfParam.PW.PPCGsbSize;


% ***********************************************************************
% Initialization
% ***********************************************************************

width = EigSol.psi.NumStateTotal();

if numEig > width
    msg = "Number of eigenvalues requested = " + num2str(numEig) + ...
          ", which is larger than the number of columns in psi = " + num2str(width);
    error(msg);
end

X = EigSol.psi;
H = EigSol.hamKS;
F = EigSol.fft;

% NOTE:
% S = ( X | W | P ) is a triplet used for PPCG.
% W is the preconditioned residual
% intermediate matrix temporarily used
% AMat = S' * H * S = S' * (AS)
% BMat = S' * S

isRestart = false;

% NOTE: never perform locking in this version


% ************************************************************************
% Main loop
% ************************************************************************

% ----------------------- Orthogonalization ---------------------------

% Option 1: Cholesky factorization 
% XTX = X' * X;
% U = chol(XTX); 
% X = X / U;

% Option 2: QR decomposition
[X, ~] = qr(X, 0);

AX = H * X;


% --------------------- Start the main loop ----------------------------
iter = 0;

% initialize size of P and AP
P = X;
AP = X;

while iter < eigMaxIter
    iter = iter + 1;
    
    if iter == 1 || isRestart
        numSet = 2;  % Steepest descent (Davidson), only use (X | W)
    else
        numSet = 3;  % Conjugate gradient, use all the triplet (X | W | P)
    end
    
    % compute the residual
    XTAX = X' * AX;
    R = AX - X * XTAX;
        
    % compute the preconditioned residual W = T * R
    W = F' * (F*R .* F.TeterPrecond);
    W = real(W);
    AW = H * W;
    
    XTW = X' * W;
    W = W - X * XTW;
    AW = AW - AX * XTW;
    
    % normalize columns of W
    normW = vecnorm(W, 2, 1);
    W = W ./ normW;
    AW = AW ./ normW;
    
    if numSet == 3
        XTP = X' * P;
        P = P - X * XTP;
        AP = AP - AX * XTP;
        
        % normalize the conjugate direction
        normP = vecnorm(P, 2, 1);
        P = P ./ normP;
        AP = AP ./ normP;
    end

    % perform the sweep
    nsb = width / sbSize;  % this should be generalized to subblocks
        
    for k = 1 : nsb
        % fetch individual columns
        idxSta = sbSize * (k - 1) + 1;
        idxEnd = idxSta + sbSize - 1;
        
        x  =  X(:, idxSta : idxEnd);
        w  =  W(:, idxSta : idxEnd);
        ax = AX(:, idxSta : idxEnd);
        aw = AW(:, idxSta : idxEnd);
        
        if numSet == 3
%             p = P(:, k : (k + sbSize - 1));
%             ap = AP(:, k : (k + sbSize -1));
            p = P(:, idxSta : idxEnd);
            ap = AP(:, idxSta : idxEnd);
        end            
        
        % compute AMat
        if numSet == 2
            xTax = x' * ax;
            xTaw = x' * aw;
            wTaw = w' * aw;
            AMat = [xTax, xTaw; 
                    xTaw', wTaw];
        end
        if numSet == 3
            xTax = x' * ax;
            xTaw = x' * aw;
            wTaw = w' * aw;
            xTap = x' * ap;
            wTap = w' * ap;
            pTap = p' * ap;
            AMat = [xTax, xTaw, xTap; 
                    xTaw', wTaw, wTap;
                    xTap', wTap', pTap];
        end
        
        % compute BMat (positive definite)
        if numSet == 2
            s = [x, w];
            BMat = s' * s;
        end
        if numSet == 3
            s = [x, w, p];
            BMat = s' * s;
        end
        
        if numSet == 3
            dim = 3 * sbSize;
        else
            dim = 2 * sbSize;
        end
        
        AMat = (AMat + AMat') / 2;
        BMat = (BMat + BMat') / 2;
        
        [V, D] = eigs(AMat(1:dim, 1:dim), BMat(1:dim, 1:dim), dim);
        
        d = real(diag(D));
        [~, id] = sort(d);  % in ascending order
        AMat = V(:, id);

%         idxSta = sbSize * (k - 1) + 1;
%         idxEnd = idxSta + sbSize - 1;
%         x  =  X(:, idxSta : idxEnd);
%         w  =  W(:, idxSta : idxEnd);
%         p  =  P(:, idxSta : idxEnd);
%         ax = AX(:, idxSta : idxEnd);
%         aw = AW(:, idxSta : idxEnd);
%         ap = AP(:, idxSta : idxEnd);
        
        cx = AMat(1 : sbSize, 1 : sbSize);
        cw = AMat((sbSize+1) : 2*sbSize, 1 : sbSize);

        if numSet == 3
            cp = AMat(2*sbSize+1 : 3*sbSize, 1 : sbSize);
            
            p = w * cw + p * cp;
            ap = aw * cw + ap * cp;
        else
            p = w * cw;
            ap = aw * cw;
        end
        
        x = x * cx + p;
        ax = ax * cx + ap;
        
        X(:, idxSta : idxEnd) = x;
        W(:, idxSta : idxEnd) = w;
        P(:, idxSta : idxEnd) = p;
        AX(:, idxSta : idxEnd) = ax;
        AW(:, idxSta : idxEnd) = aw;
        AP(:, idxSta : idxEnd) = ap;
    end        
    
    % CholeskyQR of the updated block X
    XTX = X' * X;
    U = chol(XTX); 
    X = X / U;
    
end

% ************************************************************************
% Post processing
% ************************************************************************

% Obtain the eigenvalues and eigenvectors
% XTAX should contain the matrix X' * (AX) and X is an orthonormal set

XTAX = X' * AX;

XTAX = (XTAX + XTAX') / 2;

[V, D] = eigs(XTAX, width);
d = real(diag(D));
[eigValS, id] = sort(d);
XTAX = V(:, id);

X = X * XTAX;
AX = AX * XTAX;

% compute norms of individual eigenpairs
Xtemp = AX - X .* eigValS';

resNorm = vecnorm(Xtemp, 2, 1);
resNorm = resNorm ./ max(1.0, abs(eigValS)');
resMax = max(resNorm);
resMin = min(resNorm);

% save the eigenvalues and eigenvectors back to the EigenSolverKS data
% structure

EigSol.eigVal = eigValS;
EigSol.resVal = resNorm;
EigSol.psi = X;

InfoPrint(0, "\nAfter %d iterations of PPCG \n", iter);
InfoPrint(0, "The maximum norm of the residual is %1.8e \n", resMax); 
InfoPrint(0, "The minimum norm of the residual is %1.8e \n\n", resMin); 


end