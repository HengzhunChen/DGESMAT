function scfDG = Setup(scfDG, hamDG, vecEigSol, fft, ptable)
% SCFDG/SETUP initializes SCFDG object scfDG.
%
%    scfDG = Setup(scfDG, hamDG, vecEigSol, fft, ptable) returns a SCFDG
%    object with respect to HamiltonianDG object hamdG, cell vecEigSol
%    containing EigenSolverKS object over each extended element, Fourier
%    object fft and PeriodTable object ptable.
%
%    See also SCFDG, HamiltonianDG, EigenSolverKS, Fourier, PeriodTable, 
%    ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;


% *********************************************************************
% Basic Parameters
% *********************************************************************

scfDG.domain    = esdfParam.basic.domain;
scfDG.hamDG     = hamDG;
scfDG.vecEigSol = vecEigSol;
scfDG.fft       = fft;
scfDG.ptable    = ptable;

scfDG.XCType            = esdfParam.basic.XCType;
scfDG.VDWType           = esdfParam.basic.VDWType;

scfDG.PWSolver          = esdfParam.basic.PWSolver;
scfDG.DGSolver          = esdfParam.basic.DGSolver;

scfDG.numUnusedState    = esdfParam.basic.numUnusedState;
scfDG.Tbeta             = esdfParam.basic.Tbeta;
scfDG.Tsigma            = 1.0 / scfDG.Tbeta;
scfDG.numElem           = esdfParam.DG.numElem;
scfDG.ecutWavefunction  = esdfParam.basic.ecutWavefunction;
scfDG.densityGridFactor = esdfParam.basic.densityGridFactor;
scfDG.LGLGridFactor     = esdfParam.DG.LGLGridFactor;

scfDG.periodicPotential.distancePeriodize ...
                        = esdfParam.DG.distancePeriodize;
scfDG.potentialBarrier.BarrierW = esdfParam.DG.potentialBarrierW;
scfDG.potentialBarrier.BarrierS = esdfParam.DG.potentialBarrierS;
scfDG.potentialBarrier.BarrierR = esdfParam.DG.potentialBarrierR;

% FIXME fixed ratio between the size of the extended element and the
% element
for d = 1 : dimDef()
    if scfDG.numElem(d) > 1
        scfDG.extElemRatio(d) = 3;
    else
        scfDG.extElemRatio(d) = 1;
    end
end


% ---------------------- control parameters ---------------------------

scfDG.controlVar.eigMinTolerance         = esdfParam.control.eigMinTolerance;
scfDG.controlVar.eigTolerance            = esdfParam.control.eigTolerance;
scfDG.controlVar.eigMinIter              = esdfParam.control.eigMinIter;
scfDG.controlVar.eigMaxIter              = esdfParam.control.eigMaxIter;
scfDG.controlVar.scfInnerTolerance       = esdfParam.control.scfInnerTolerance;
scfDG.controlVar.scfInnerMinIter         = esdfParam.control.scfInnerMinIter;
scfDG.controlVar.scfInnerMaxIter         = esdfParam.control.scfInnerMaxIter;
scfDG.controlVar.scfOuterTolerance       = esdfParam.control.scfOuterTolerance;
scfDG.controlVar.scfOuterMinIter         = esdfParam.control.scfOuterMinIter;
scfDG.controlVar.scfOuterMaxIter         = esdfParam.control.scfOuterMaxIter;
scfDG.controlVar.scfOuterEnergyTolerance = esdfParam.control.scfOuterEnergyTolerance;
scfDG.controlVar.SVDBasisTolerance       = esdfParam.control.SVDBasisTolerance;
scfDG.controlVar.isPWeigToleranceDynamic = esdfParam.userOption.general.isPWeigTolDynamic;


% ------------------- Chebyshev Filter --------------------------------

% Chebyshev filtering related parameters for PWDFT on extended element
scfDG.CheFSIPW = esdfParam.PW.CheFSI;

% Chebyshev filerting related parameters for DG solver
if scfDG.DGSolver == "CheFSI"

    scfDG.CheFSIDG.firstFilterOrder   = esdfParam.DG.CheFSI.firstFilterOrder;
    scfDG.CheFSIDG.firstCycleNum      = esdfParam.DG.CheFSI.firstCycleNum;
    scfDG.CheFSIDG.secondOuterIter    = esdfParam.DG.CheFSI.secondOuterIter;
    scfDG.CheFSIDG.secondFilterOrder  = esdfParam.DG.CheFSI.secondFilterOrder;
    scfDG.CheFSIDG.secondCycleNum     = esdfParam.DG.CheFSI.secondCycleNum;
    scfDG.CheFSIDG.generalFilterOrder = esdfParam.DG.CheFSI.generalFilterOrder;
    scfDG.CheFSIDG.generalCycleNum    = esdfParam.DG.CheFSI.generalCycleNum;
    
    scfDG.ionDyn.isChebyInIonDyn = 0;
    scfDG.ionDyn.ionDynIter = 0;

    % ChebyShev polynomial filtered complementary subspace iteration
    % Only accessed if CheFSI is in use
    scfDG.CheFSIDG.isUseCompSubspace = esdfParam.DG.CheFSI.isUseCompSubspace;  % default 0

    % Safeguard to ensure that CS strategy is called only after at least 
    % one general CheFSI cycle has been called
    % This allows the initial guess vectors to be copied
    if scfDG.CheFSIDG.isUseCompSubspace && scfDG.CheFSIDG.secondOuterIter < 2
        scfDG.CheFSIDG.secondOuterIter = 2;
    end

    CompSubspaceParam = esdfParam.DG.CheFSI.CompSubpace;
    scfDG.CheFSIDG.CompSubspace.nStates = ...
        CompSubspaceParam.nStates;  % Defaults to a fraction of extra states
    scfDG.CheFSIDG.CompSubspace.ionIterRegularChebyFreq = ...
        CompSubspaceParam.ionIterRegularChebyFreq;  % Defaults to 20
    scfDG.CheFSIDG.CompSubspace.biggerGridDimFac = ...
        CompSubspaceParam.biggerGridDimFac;  % Defaults to 1;

    % LOBPCG for top states option
    scfDG.CheFSIDG.CompSubspace.lobpcgIter = ...
        CompSubspaceParam.lobpcgIter;  % Default = 15
    scfDG.CheFSIDG.CompSubspace.lobpcgTol = ...
        CompSubspaceParam.lobpcgTol;  % Default = 1e-8

    % CheFSI for top states option
    scfDG.CheFSIDG.CompSubspace.isHMatTopStatesUseCheby = ...
        CompSubspaceParam.isHMatTopStatesUseCheby;
    scfDG.CheFSIDG.CompSubspace.HMatFilterOrder = ...
        CompSubspaceParam.HMatFilterOrder; 
    scfDG.CheFSIDG.CompSubspace.HMatCycleNum = ...
        CompSubspaceParam.HMatCycleNum; 
    scfDG.CheFSIDG.CompSubspace.HMatdeltaFudge = 0.0;

    scfDG.CheFSIDG.CompSubspace.nSolve = ...
        hamDG.NumExtraState() + scfDG.CheFSIDG.CompSubspace.nStates;     
else
    scfDG.CheFSIDG.isUseCompSubspace = 0;
end

% -------------- end of ChebyShev polynomial filter -------------------


% -----------------------  Mixing  -------------------------------------

scfDG.mix.mixMaxDim     = esdfParam.basic.mixMaxDim;
scfDG.mix.mixVariable   = esdfParam.basic.mixVariable;
scfDG.mix.mixType       = esdfParam.basic.mixType;
scfDG.mix.mixStepLength = esdfParam.basic.mixStepLength;

numElem = scfDG.numElem;
scfDG.mix.mixOuterSave = cell(numElem);
scfDG.mix.mixInnerSave = cell(numElem);

scfDG.mix.dfOuterMat = cell(numElem);
scfDG.mix.dvOuterMat = cell(numElem);
scfDG.mix.dfInnerMat = cell(numElem);
scfDG.mix.dvInnerMat = cell(numElem);

scfDG.vtotLGLSave = cell(numElem);

mixMaxDim = scfDG.mix.mixMaxDim;
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            numUniformGridTotalFine = prod( hamDG.grid.numUniformGridElemFine );
            
            scfDG.mix.mixOuterSave{i, j, k} = zeros(numUniformGridTotalFine, 1);
            scfDG.mix.mixInnerSave{i, j, k} = zeros(numUniformGridTotalFine, 1);
            
            scfDG.mix.dfOuterMat{i, j, k} = zeros(numUniformGridTotalFine, mixMaxDim);
            scfDG.mix.dvOuterMat{i, j, k} = zeros(numUniformGridTotalFine, mixMaxDim);
            scfDG.mix.dfInnerMat{i, j, k} = zeros(numUniformGridTotalFine, mixMaxDim);
            scfDG.mix.dvInnerMat{i, j, k} = zeros(numUniformGridTotalFine, mixMaxDim);
            
            numLGLGridTotal = prod( hamDG.grid.numLGLGridElem );
            scfDG.vtotLGLSave{i, j, k} = zeros(numLGLGridTotal, 1);
        end
    end
end


% ------------------------  Smearing  -----------------------------------

% choice of smearing scheme : Fermi-Dirac (FD) or Gaussian-Broadening (GB)
% or Methfessel-Paxton (MP) 
scfDG.smearing.SmearingScheme = esdfParam.basic.smearingScheme;
if scfDG.smearing.SmearingScheme == "GB"
    scfDG.smearing.MPsmearingOrder = 0;
elseif scfDG.smearing.SmearingScheme == "MP"
    scfDG.smearing.MPsmearingOrder = 2;
else
    scfDG.smearing.MPsmearingOrder = -1;  % for safety
end


% ------------------------  Restart  ------------------------------------

% restart the density in the global domain
scfDG.restart.DensityFileName = "DENSITY.mat";
% restart the wavefunctions in the extended element
scfDG.restart.WfnExtElemFileName = "WFNEXTELEM.mat";


% --------------  Ionic iteration related parameters  -------------------

% Ionic iteration number
scfDG.ionDyn.IonDynIter = 0;  
% Whether to use energy based SCF convergence
scfDG.ionDyn.isUseEnergySCFconvergence = 0;  
 % Tolerance for SCF total energy for energy based SCF convergence
scfDG.ionDyn.MDscfEtotDiffTol = esdfParam.ionDyn.MDscfEtotDiffTol; 
% Tolerance for SCF band energy for energy based SCF convergence
scfDG.ionDyn.MDscfEbandDiffTol = esdfParam.ionDyn.MDscfEbandDiffTol;  

scfDG.ionDyn.MDscfEtot = 0.0;
scfDG.ionDyn.MDscfEtotOld = 0.0;
scfDG.ionDyn.MDscfEtotDiff = 0.0;
scfDG.ionDyn.MDscfEband = 0.0;
scfDG.ionDyn.MDscfEbandOld = 0.0; 
scfDG.ionDyn.MDscfEbandDiff = 0.0;


% --------------  end of basic SCFDG paramters  -----------------------


% *********************************************************************
% Initialization
% *********************************************************************

numElem = scfDG.numElem;

% assume the initial error is O(1)
scfDG.scfOuterNorm = 1.0;
scfDG.scfInnerNorm = 1.0;

% initial value
scfDG.efreeDifPerAtom = 100.0;


% -------------------------- density ---------------------------------

if esdfParam.userOption.general.isRestartDensity
    reStartData = load(scfDG.restart.DensityFileName);
    density = reStartData.density;
    sumDensity = 0;
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                if size(density{i, j, k}) ~= size(scfDG.hamDG.density{i, j, k})
                    error('The size of restarting density does not match with current setup');
                end
                sumDensity = sumDensity + sum(density{i, j, k});
            end
        end
    end
    scfDG.hamDG.density = density;

    InfoPrint(0, "Restart density. Sum of density      = ", ...
              sumDensity * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine() );
    
else  % using zero intial guess
    if esdfParam.userOption.general.isUseAtomDensity
        atomDensity = scfDG.hamDG.CalculateAtomDensity(scfDG.ptable);        
        scfDG.hamDG.atomDensity = atomDensity;
        scfDG.hamDG.density = atomDensity;
        
    else  % use pseudoCharges
        % Initialize the electron density using the pseudocharge, make sure
        % the pseudocharge is initialized
        sumDensity = 0;
        sumPseudoCharge = 0;
        EPS = 1e-6;
        
        % make sure that the electron density is positive
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    pseudoCharge = scfDG.hamDG.pseudoCharge{i, j, k};
                    idxnz = pseudoCharge > EPS;
                    density = pseudoCharge;
                    density(~idxnz) = EPS;
                    scfDG.hamDG.density{i, j, k} = density;
                    sumDensity = sumDensity + sum(density);
                    sumPseudoCharge = sumPseudoCharge + sum(pseudoCharge);
                end
            end
        end
        
        % rescale the density
        InfoPrint( 0, "Initial density. Sum of density      = ",  ...
            sumDensity .* scfDG.domain.Volume() ./ scfDG.domain.NumGridTotalFine() );
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    scfDG.hamDG.density{i, j, k} = scfDG.hamDG.density{i, j, k} .* ...
                        (sumPseudoCharge / sumDensity);
                end
            end
        end
    end
end     

% ----------------------- end of density -------------------------------


% -------------- wavefunctions in extended element ---------------------

if esdfParam.userOption.general.isRestartWfn
    reStartData = load(scfDG.restart.WfnExtElemFileName);
    wfnExtElem = reStartData.wfnExtElem;
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                if size(wfnExtElem{i, j, k}) ~= size(scfDG.vecEigSol{i, j, k}.psi.wavefun)
                    error('The size of restarting wavefun does not match with current setup');
                end
                scfDG.vecEigSol{i, j, k}.psi.wavefun = wfnExtElem{i, j, k};    
            end
        end
    end
else
    % use random initial guess for basis functions in the extended element
    InfoPrint(0, "Initial basis functions with random guess.");
end

% -------------------- end of wavefunctions ----------------------------


% ------------ transfer matrix from uniform grid to LGL grid -----------

% Generate the transfer matrix from the periodic uniform grid on each
% extended element to LGL grid. 
% The interpolation must be performed through a fine Fourier grid 
% (uniform grid) and then interpolate to the LGL grid.

scfDG.PeriodicUniformToLGLMat = cell(dimDef(), 1);
scfDG.PeriodicUniformFineToLGLMat = cell(dimDef(), 1);
scfDG.PeriodicGridExtElemToGridElemMat = cell(dimDef(), 1);

dmExtElem = scfDG.vecEigSol{1, 1, 1}.fft.domain;
dmElem = Domain();
for d = 1 : dimDef()
    dmElem.length(d) = scfDG.domain.length(d) / numElem(d);
    dmElem.numGrid(d) = scfDG.domain.numGrid(d) / numElem(d);
    dmElem.numGridFine(d) = scfDG.domain.numGridFine(d) / numElem(d);
    % posStart relative to the extended element
    dmExtElem.posStart(d) = 0;
    if numElem(d) > 1
        dmElem.posStart(d) = dmElem.length(d);
    else
        dmElem.posStart(d) = 0;
    end
end

numLGL             = hamDG.grid.numLGLGridElem;
numUniform         = dmExtElem.numGrid;
numUniformFine     = dmExtElem.numGridFine;
numUniformFineElem = dmElem.numGridFine;
lengthUniform      = dmExtElem.length;
    
LGLGrid             = LGLMesh(dmElem, numLGL);
UniformGrid         = UniformMesh(dmExtElem);
UniformGridFine     = UniformMeshFine(dmExtElem);
UniformGridFineElem = UniformMeshFine(dmElem);

for d = 1 : dimDef()    
    ng = numUniform(d);
    grid = (0:ng-1) - ng *( (0:ng-1) > ng/2 );
    KGrid = grid * 2 * pi ./ lengthUniform(d);  % row vector
    
    tempMat = LGLGrid{d}' - UniformGrid{d};
    tempTns = repmat(tempMat, [1, 1, numUniform(d)]); 
    % size [numLGL(d), numUniform(d), numUniform(d)]
    tempKGrid = repmat(KGrid, [numLGL(d), 1, numUniform(d)]);
    tempKGrid = permute(tempKGrid, [1, 3, 2]);  
    % size [numLGL(d), numUniform(d), numUniform(d)], i-th page is all KGrid(i)
    localMat = sum( cos(tempTns .* tempKGrid), 3 ) ./ numUniform(d);  
    % size [numLGL(d), numUniform(d)]

    tempMatFineElem = UniformGridFineElem{d}' - UniformGrid{d};
    tempTnsFineElem = repmat(tempMatFineElem, [1, 1, numUniform(d)]);
    tempKGrid = repmat(KGrid, [numUniformFineElem(d), 1, numUniform(d)]);
    tempKGrid = permute(tempKGrid, [1, 3, 2]);
    localMatFineElem = sum( cos(tempTnsFineElem .* tempKGrid), 3 ) ./ numUniform(d);
    % size [numUniformFineElem(d), numUniform(d)]
    
    scfDG.PeriodicUniformToLGLMat{d} = localMat;
    scfDG.PeriodicGridExtElemToGridElemMat{d} = localMatFineElem;
    
    % used for reference
%     for i = 1 : numLGL(d)
%         for j = 1 : numUniform(d)
%             localMat(i, j) = 0.0;
%             for k = 1 : numUniform(d)
%                 localMat(i, j) = localMat(i, j) + ...
%                     cos( KGrid(k) * ( LGLGrid(d)(i) - UniformGrid(d)(j) ) ) / numUniform(d);
%             end
%         end
%     end     
end

for d = 1 : dimDef()
    ng = numUniformFine(d);
    grid = (0:ng-1) - ng *( (0:ng-1) > ng/2 );
    KGridFine = grid * 2 * pi ./ lengthUniform(d);  % row vector

    tempMatFine = LGLGrid{d}' - UniformGridFine{d};
    tempTnsFine = repmat(tempMatFine, [1, 1, numUniformFine(d)]); 
    % size [numLGL(d), numUniformFine(d), numUniformFine(d)]
    tempKGridFine = repmat(KGridFine, [numLGL(d), 1, numUniformFine(d)]);
    tempKGridFine = permute(tempKGridFine, [1, 3, 2]);  
    % size [numLGL(d), numUniformFine(d), numUniformFine(d)], i-th page is all KGridFine(i)    
    localMatFine = sum( cos(tempTnsFine .* tempKGridFine), 3 ) ./ numUniformFine(d);  
    % size [numLGL(d), numUniformFine(d)]
    
    scfDG.PeriodicUniformFineToLGLMat{d} = localMatFine;
end

% ----------------- end of grid transfer matrix -------------------------


% ---------------------- potential barrier ---------------------------

% Whether to apply potential barrier in the extended element. CANNOT be
% used together with periodization option.
if esdfParam.userOption.DG.isPotentialBarrier
    DIM = dimDef();
    scfDG.potentialBarrier.vBarrier = cell(DIM, 1);
    barrierS = scfDG.potentialBarrier.BarrierS;
    barrierW = scfDG.potentialBarrier.BarrierW;
    barrierR = scfDG.potentialBarrier.BarrierR;

    dmExtElem = scfDG.vecEigSol{1, 1, 1}.fft.domain;
    gridpos = UniformMeshFine(dmExtElem);

    for d = 1 : DIM
        dmlength = dmExtElem.length(d);
        numGridFine = dmExtElem.numGridFine(d);
        posStart = dmExtElem.posStart(d);
        center = posStart + dmlength / 2;

        % FIXME
        EPS = 1.0;  % for stability reason
        vBarrier = zeros(numGridFine, 1);

        dist = abs(gridpos{d} - center);
        % only apply the barrier for region outside barrierR
        idx = dist > barrierR;
        vBarrier(idx) = barrierS ...
            .* exp( -barrierW ./ (dist(idx) - barrierR) ) ...
            ./ (dist(idx) - dmlength/2 - EPS).^2 ; 
        scfDG.potentialBarrier.vBarrier{d} = vBarrier;
    end
end
            
% ------------------ end of potential barrier -------------------------


% --------------------- periodize the potential ------------------------

% Whether to periodize the potential in the extended element. CANNOT be
% used together with barrier option.
if esdfParam.userOption.DG.isPeriodizePotential
    DIM = dimDef();
    scfDG.periodicPotential.vBubble = cell(DIM, 1);
    distancePeriodize = scfDG.periodicPotential.distancePeriodize;
    
    dmExtElem = scfDG.vecEigSol{i, j, k}.fft.domain;
    gridpos = UniformMeshFine( dmExtElem );

    for d = 1 : DIM
        dmlength = dmExtElem.length(d);
        numGridFine = dmExtElem.numGridFine(d);
        posStart = dmExtElem.posStart(d);

        % FIXME
        EPS = 0.2;  % criterion for distancePeriodize
        vBubble = ones(numGridFine, 1);

        if distancePeriodize(d) > EPS
            lb = posStart + distancePeriodize(d);
            rb = posStart + dmlength - distancePeriodize(d);
            idxr = gridpos{d} > rb;
            vBubble(idxr) = Smoother( (gridpos{d}(idxr) - rb) ...
                ./ (distancePeriodize(d) - EPS) ); 
            idxl = gridpos{d} < lb;
            vBubble(idxl) = Smoother( (lb - gridpos{d}(idxl)) ...
                ./ (distancePeriodize(d) - EPS) );
        end

        scfDG.periodicPotential.vBubble{d} = vBubble;
    end
end

% ------------------ end of periodize the potential ---------------------


end