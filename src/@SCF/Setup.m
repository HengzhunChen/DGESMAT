function scf = Setup(scf, eigSol)
% SCF/SETUP initializes SCF object scf with EigenSolverKS object eigSol and
%    global variable esdfParam.
%
%    See also SCF, EigenSolverKS, ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

% --------------------  basic parameters  -------------------------------

scf.Tbeta  = esdfParam.basic.Tbeta;
scf.eigSol = eigSol;

scf.PWSolver = esdfParam.basic.PWSolver;
% Chebyshev Filtering related parameters
scf.CheFSI = esdfParam.PW.CheFSI;
scf.isChebyInIonDyn = 0;


% ---------------------  control variables ------------------------------

% Note: for PW SCF there is no inner loop. Use the parameter value
% for the outer SCF loop only.
scf.controlVar.eigTolerance          = esdfParam.control.eigTolerance;
scf.controlVar.eigMinTolerance       = esdfParam.control.eigMinTolerance;
scf.controlVar.eigMaxIter            = esdfParam.control.eigMaxIter;
scf.controlVar.scfTolerance          = esdfParam.control.scfOuterTolerance;
scf.controlVar.scfMaxIter            = esdfParam.control.scfOuterMaxIter;
scf.controlVar.scfPhiMaxIter         = esdfParam.hybrid.scfPhiMaxIter;
scf.controlVar.scfPhiTolerance       = esdfParam.hybrid.scfPhiTolerance;
scf.controlVar.isEigToleranceDynamic = esdfParam.userOption.general.isPWeigTolDynamic;


% -------------------------  mixing  ----------------------------------

scf.mixing.mixMaxDim     = esdfParam.basic.mixMaxDim;
scf.mixing.mixType       = esdfParam.basic.mixType;
scf.mixing.mixStepLength = esdfParam.basic.mixStepLength;

ntotFine = esdfParam.basic.domain.NumGridTotalFine();
scf.vtotNew = zeros(ntotFine, 1);
scf.mixing.dfMat = zeros(ntotFine, scf.mixing.mixMaxDim);
scf.mixing.dvMat = zeros(ntotFine, scf.mixing.mixMaxDim);


% -----------------------  restart  ----------------------------------

scf.restart.DensityFileName   = 'DENSITY.mat';
scf.restart.WfnFileName       = 'WAVEFUN.mat';


% -------------------------- Density ----------------------------------

if esdfParam.userOption.general.isRestartDensity
    reStartData = load(scf.restart.DensityFileName);
    scf.eigSol.hamKS.density = reStartData.density;
    InfoPrint(0, 'Density restarted from file %s \n', scf.restart.DensityFileName);

else % using the zero initial guess
    if esdfParam.userOption.general.isUseAtomDensity
        % use the superposition of atomic density as the initial guess for density 
        InfoPrint(0, 'Use superposition of atomic density as initial guess for electron density. \n');
        
        ptable = eigSol.hamKS.ptable;
        atomDensity = scf.eigSol.hamKS.CalculateAtomDensity(ptable);                
        scf.eigSol.hamKS.density = atomDensity;
        
    else
        % Start from pseudocharge, usually this is not a very good idea
        % make sure the pseudocharge is initialized
        InfoPrint(0, 'Generating initial density through linear combination of pseudocharges.\n');
        
        domain = esdfParam.basic.domain;
        pseudoCharge = scf.eigSol.hamKS.pseudoCharge;
        EPS = 1e-6;
        
        % make sure that the electron density is positive
        density = max(pseudoCharge, EPS);
        sum0 = sum(density);
        sum1 = sum(pseudoCharge);
        
        InfoPrint( 0, "Initial density. Sum of density       = ", ...
            sum0 * domain.Volume() / domain.NumGridTotalFine() );
        
        % rescale the density
        density = density .* sum1 ./ sum0;
        scf.eigSol.hamKS.density = density;
        
        InfoPrint( 0, "Rescaled density. Sum of density      = ", ... 
            sum1 * domain.Volume() / domain.NumGridTotalFine() );
    end    
end


% -------------------------- wave function ----------------------------

if ~esdfParam.userOption.general.isRestartWfn
    % randomized input from outside
    % setup the occupation rate by aufbau principle 
    % (needed for hybrid functional calculation)
    npsi = eigSol.psi.NumStateTotal();
    occ = ones(npsi, 1);
    scf.eigSol.hamKS.occupationRate = occ;
else
    reStartData = load(scf.restart.WfnFileName);

    scf.eigSol.psi.wavefun = reStartData.wavefun;
    scf.eigSol.hamKS.occupationRate = reStartData.occupationRate;

    InfoPrint(0, 'Wavefunction restarted from file %s \n', scf.restart.WfnFileName);
end
    

% -------------------------- XC functional ----------------------------

scf.controlVar.isCalculateGradRho = false;
if contains(esdfParam.basic.XCType, "GGA")
    scf.controlVar.isCalculateGradRho = true;
end


% external force when needed


end