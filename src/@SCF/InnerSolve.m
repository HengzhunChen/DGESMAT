function scf = InnerSolve(scf, iter)
% SCF/INNERSOLVE inner process of the iter-th SCF iteration, including 
%    solve the eigenvalue problem and update energy and potential.
%
%    See also SCF, EigenSolverKS, HamiltonianKS.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


controlVar = scf.controlVar; 

% **********************************************************************
%                   solve the eigenvalue problem                               
% **********************************************************************

if controlVar.isEigToleranceDynamic
    % Dynamic strategy to control the tolerance
    if iter == 1
        eigTolNow = 1e-2;
    else
        eigTolNow = max( min(scf.scfNorm*1e-2, 1e-2), controlVar.eigTolerance );
    end
else
    % Static strategy to control the tolerance
    eigTolNow = controlVar.eigTolerance;
end

eigMaxIter = controlVar.eigMaxIter;
eigMinTolerance = controlVar.eigMinTolerance;

numEig = scf.eigSol.psi.NumStateTotal();

% if not use Chebyshev
if scf.PWSolver ~= "CheFSI"
    InfoPrint(0, 'The current tolerance used by the eigensolver is %1.8e \n', eigTolNow);
    InfoPrint(0, 'The target number of converged eigenvectors is %d \n', numEig);
end

eigSol = scf.eigSol;

if scf.PWSolver == "LOBPCG"
    % Use LOBPCG
    eigSol = LOBPCGSolve(eigSol, numEig, eigMaxIter, eigMinTolerance, eigTolNow);
elseif scf.PWSolver == "PPCG"
    % Use PPCG
    eigSol = PPCGSolve(eigSol, numEig, eigMaxIter, scf.PPCGsbSize);
elseif scf.PWSolver == "eigs"
    % Use eigs() function in MATLAB
    eigSol = eigsSolve(eigSol, numEig, eigMaxIter, eigTolNow);
elseif scf.PWSolver == "CheFSI"
    % Use CheFSI
    InfoPrint(0, ' CheFSI in PWDFT working on static schedule. \n');
    % use CheFSI or LOBPCG on first step
    if iter <= 1
        if scf.CheFSI.firstCycleNum <= 0
            eigSol = LOBPCGSolve(eigSol, numEig, ...
                eigMaxIter, eigMinTolerance, eigTolNow);
        else
            eigSol = ChebyStepFirst(eigSol, numEig, ...
                scf.CheFSI.firstCycleNum, scf.CheFSI.firstFilterOrder);
        end
    else
        eigSol = ChebyStepGeneral(eigSol, numEig, scf.CheFSI.generalFilterOrder);
    end
else
    error('Not supported PWSolver type.');
end

scf.eigSol = eigSol;
scf.eigSol.hamKS.eigVal = scf.eigSol.eigVal;


% **********************************************************************
%                 compute energy and potential, etc                              
% **********************************************************************

[scf.eigSol.hamKS.occupationRate, scf.fermi] = ...
        scf.CalculateOccupationRate(scf.eigSol.eigVal);

% Calculate the Harris energy before updating the density
scf.EfreeHarris = scf.CalculateHarrisEnergy();

hamKS = scf.eigSol.hamKS;

[scf.totalCharge, hamKS.density] = CalculateDensity(...
            hamKS, scf.eigSol.psi, hamKS.occupationRate);
if controlVar.isCalculateGradRho
    hamKS.gradDensity = CalculateGradDensity(hamKS);
end

[scf.Exc, hamKS.epsxc, hamKS.vxc] = CalculateXC(hamKS);
hamKS.vhart = CalculateHartree(hamKS);
% NO external potential
scf.vtotNew = CalculateVtot(hamKS);

scf.eigSol.hamKS = hamKS;


end