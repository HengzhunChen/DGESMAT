function scf = Setup(scf, esdfParam, eigSol)
% SCF/SETUP initializes SCF object scf with EigenSolverKS object eigSol and
%    ESDFInputParam object esdfParam.
%
%    See also SCF, EigenSolverKS, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% --------------------- user options ----------------------------------

userOptionPW = esdfParam.userOption.PW;

scf.userOption.isOutputWavefun    = userOptionPW.isOutputWavefun;
scf.userOption.isOutputDensity    = userOptionPW.isOutputDensity;
scf.userOption.isOutputPotential  = userOptionPW.isOutputPotential;
scf.userOption.isOutputAtomStruct = userOptionPW.isOutputAtomStruct;

scf.userOption.isUseVLocal = esdfParam.userOption.general.isUseVLocal;


% --------------- IO data file names -------------------------------

dataFileIO = esdfParam.dataFileIO;

scf.dataFileIO.wavefunPW    = dataFileIO.wavefunPW;
scf.dataFileIO.densityPW    = dataFileIO.densityPW;
scf.dataFileIO.potentialPW  = dataFileIO.potentialPW;
scf.dataFileIO.atomStructPW = dataFileIO.atomStructPW;

scf.dataFileIO.restartDensity = dataFileIO.restartDensity;
scf.dataFileIO.restartWfn     = dataFileIO.restartWfn;


% --------------------  basic parameters  -------------------------------

scf.Tbeta  = esdfParam.basic.Tbeta;
scf.eigSol = eigSol;

scf.PWSolver = esdfParam.basic.PWSolver;
scf.PPCGsbSize = esdfParam.PW.PPCGsbSize;
scf.CheFSI = esdfParam.PW.CheFSI;  % Chebyshev Filtering parameters


% ---------------------  control variables ------------------------------

% Note: for PW SCF there is no inner loop. Use the parameter value
% for the outer SCF loop only.
scf.controlVar.eigTolerance          = esdfParam.control.eigTolerance;
scf.controlVar.eigMinTolerance       = esdfParam.control.eigMinTolerance;
scf.controlVar.eigMaxIter            = esdfParam.control.eigMaxIter;
scf.controlVar.scfTolerance          = esdfParam.control.scfOuterTolerance;
scf.controlVar.scfMaxIter            = esdfParam.control.scfOuterMaxIter;
scf.controlVar.isEigToleranceDynamic = esdfParam.userOption.general.isPWeigTolDynamic;


% -------------------------  mixing  ----------------------------------

scf.mixing.mixMaxDim     = esdfParam.basic.mixMaxDim;
scf.mixing.mixType       = esdfParam.basic.mixType;
scf.mixing.mixStepLength = esdfParam.basic.mixStepLength;

ntotFine = esdfParam.basic.domain.NumGridTotalFine();
scf.vtotNew = zeros(ntotFine, 1);
scf.mixing.dfMat = zeros(ntotFine, scf.mixing.mixMaxDim);
scf.mixing.dvMat = zeros(ntotFine, scf.mixing.mixMaxDim);


% -------------------------- Density ----------------------------------

if esdfParam.userOption.general.isRestartDensity
    reStartData = load(scf.dataFileIO.restartDensity);
    scf.eigSol.hamKS.density = reStartData.density;
    InfoPrint(0, 'Density restarted from file %s \n', scf.dataFileIO.restartDensity);

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
    reStartData = load(scf.dataFileIO.restartWfn);

    scf.eigSol.psi.wavefun = reStartData.wavefun;
    scf.eigSol.hamKS.occupationRate = reStartData.occupationRate;

    InfoPrint(0, 'Wavefunction restarted from file %s \n', scf.dataFileIO.restartWfn);
end
    

% -------------------------- XC functional ----------------------------

scf.controlVar.isCalculateGradRho = false;
if contains(esdfParam.basic.XCType, "GGA")
    scf.controlVar.isCalculateGradRho = true;
end


% external force when needed


end