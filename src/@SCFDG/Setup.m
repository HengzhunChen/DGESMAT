function scfDG = Setup(scfDG, esdfParam, hamDG, vecEigSol, fft, ptable)
% SCFDG/SETUP initializes SCFDG object scfDG.
%
%    scfDG = Setup(scfDG, esdfParam, hamDG, vecEigSol, fft, ptable) 
%    returns a SCFDG object with respect to ESDFInputParam object esdfParam, 
%    HamiltonianDG object hamdG, cell vecEigSol containing EigenSolverKS 
%    object over each extended element, Fourier object fft and PeriodTable 
%    object ptable.
%
%    See also SCFDG, HamiltonianDG, EigenSolverKS, Fourier, PeriodTable, 
%    ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% --------------------------- user options ----------------------------------

userOptionDG = esdfParam.userOption.DG;

scfDG.userOption.isOutputDensity        = userOptionDG.isOutputDensity;
scfDG.userOption.isOutputPotential      = userOptionDG.isOutputPotential;
scfDG.userOption.isOutputAtomStruct     = userOptionDG.isOutputAtomStruct;
scfDG.userOption.isOutputALBElemUniform = userOptionDG.isOutputALBElemUniform;
scfDG.userOption.isOutputALBElemLGL     = userOptionDG.isOutputALBElemLGL;
scfDG.userOption.isOutputWfnExtElem     = userOptionDG.isOutputWfnExtElem;
scfDG.userOption.isOutputPotExtElem     = userOptionDG.isOutputPotExtElem;

scfDG.userOption.isPeriodizePotential = userOptionDG.isPeriodizePotential;
scfDG.userOption.isCalculateAPosterioriEachSCF = userOptionDG.isCalculateAPosterioriEachSCF;


% --------------------------- IO data file names -----------------------------

dataFileIO = esdfParam.dataFileIO;

scfDG.dataFileIO.densityDG    = dataFileIO.densityDG;
scfDG.dataFileIO.potentialDG  = dataFileIO.potentialDG;
scfDG.dataFileIO.atomStructDG = dataFileIO.atomStructDG;

scfDG.dataFileIO.albElemUniform = dataFileIO.albElemUniform;
scfDG.dataFileIO.albElemLGL     = dataFileIO.albElemLGL;

scfDG.dataFileIO.wfnExtElem = dataFileIO.wfnExtElem;
scfDG.dataFileIO.potExtElem = dataFileIO.potExtElem;

scfDG.dataFileIO.restartWfn     = dataFileIO.restartWfn;
scfDG.dataFileIO.restartDensity = dataFileIO.restartDensity;



% *********************************************************************
% Basic Parameters
% *********************************************************************

scfDG.domain    = esdfParam.basic.domain;
scfDG.hamDG     = hamDG;
scfDG.vecEigSol = vecEigSol;
scfDG.fft       = fft;
scfDG.ptable    = ptable;
scfDG.bufferSize = esdfParam.DG.bufferSize;

scfDG.XCType            = esdfParam.basic.XCType;
scfDG.VDWType           = esdfParam.basic.VDWType;

scfDG.PWSolver          = esdfParam.basic.PWSolver;
scfDG.DGSolver          = esdfParam.basic.DGSolver;

scfDG.numUnusedState    = esdfParam.basic.numUnusedState;
scfDG.Tbeta             = esdfParam.basic.Tbeta;
scfDG.Tsigma            = 1.0 / scfDG.Tbeta;
scfDG.numElem           = esdfParam.DG.numElem;

scfDG.periodicPotential.distancePeriodize ...
                        = esdfParam.DG.distancePeriodize;


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

% Chebyshev filtering parameters for PW subproblem on extended element
scfDG.CheFSIPW = esdfParam.PW.CheFSI;
% PPCG parameter for PW subproblem on extended element
scfDG.PPCGsbSize = esdfParam.PW.PPCGsbSize;

% -----------------------  Mixing  -------------------------------------

scfDG.mix.mixMaxDim     = esdfParam.basic.mixMaxDim;
scfDG.mix.mixVariable   = esdfParam.basic.mixVariable;
scfDG.mix.mixType       = esdfParam.basic.mixType;
scfDG.mix.mixStepLength = esdfParam.basic.mixStepLength;

numElem = scfDG.numElem;
numElemTotal = prod(numElem);
scfDG.mix.mixOuterSave = cell(numElemTotal, 1);
scfDG.mix.mixInnerSave = cell(numElemTotal, 1);

scfDG.mix.dfOuterMat = cell(numElemTotal, 1);
scfDG.mix.dvOuterMat = cell(numElemTotal, 1);
scfDG.mix.dfInnerMat = cell(numElemTotal, 1);
scfDG.mix.dvInnerMat = cell(numElemTotal, 1);

scfDG.vtotLGLSave = cell(numElemTotal, 1);

mixMaxDim = scfDG.mix.mixMaxDim;

for elemIdx = 1 : numElemTotal
    numUniformGridTotalFine = prod( hamDG.grid.numUniformGridElemFine );
    
    scfDG.mix.mixOuterSave{elemIdx} = zeros(numUniformGridTotalFine, 1);
    scfDG.mix.mixInnerSave{elemIdx} = zeros(numUniformGridTotalFine, 1);
    
    scfDG.mix.dfOuterMat{elemIdx} = zeros(numUniformGridTotalFine, mixMaxDim);
    scfDG.mix.dvOuterMat{elemIdx} = zeros(numUniformGridTotalFine, mixMaxDim);
    scfDG.mix.dfInnerMat{elemIdx} = zeros(numUniformGridTotalFine, mixMaxDim);
    scfDG.mix.dvInnerMat{elemIdx} = zeros(numUniformGridTotalFine, mixMaxDim);
    
    numLGLGridTotal = prod( hamDG.grid.numLGLGridElem );
    scfDG.vtotLGLSave{elemIdx} = zeros(numLGLGridTotal, 1);
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


% --------------  end of basic SCFDG paramters  -----------------------


% *********************************************************************
% Initialization
% *********************************************************************

numElem = scfDG.numElem;
numElemTotal = prod(numElem);

% assume the initial error is O(1)
scfDG.scfOuterNorm = 1.0;
scfDG.scfInnerNorm = 1.0;

% initial value
scfDG.efreeDifPerAtom = 100.0;


% -------------------------- density ---------------------------------

if esdfParam.userOption.general.isRestartDensity
    density = cell(numElemTotal);
    sumDensity = 0;
    restartFilePrefix = scfDG.dataFileIO.restartDensity;
    for elemIdx = 1 : numElemTotal
        fileName = restartFilePrefix + "_" + num2str(elemIdx) + ".mat";
        restartData = load(fileName);
        densityElem = restartData.densityElem;
        if size(densityElem) ~= size(scfDG.hamDG.density{elemIdx})
            error('The size of restarting density does not match with current setup');
        end
        density{elemIdx} = densityElem;
        sumDensity = sumDensity + sum(density{elemIdx});
    end
    scfDG.hamDG.density = density;

    InfoPrint(0, "Restart density. Sum of density      = ", ...
              sumDensity * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine() );
    
else  % using zero initial guess
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
        for elemIdx = 1 : numElemTotal
            pseudoCharge = scfDG.hamDG.pseudoCharge{elemIdx};
            idxnz = pseudoCharge > EPS;
            density = pseudoCharge;
            density(~idxnz) = EPS;
            scfDG.hamDG.density{elemIdx} = density;
            sumDensity = sumDensity + sum(density);
            sumPseudoCharge = sumPseudoCharge + sum(pseudoCharge);
        end
        
        % rescale the density
        InfoPrint( 0, "Initial density. Sum of density      = ",  ...
            sumDensity .* scfDG.domain.Volume() ./ scfDG.domain.NumGridTotalFine() );
        for elemIdx = 1 : numElemTotal
            scfDG.hamDG.density{elemIdx} = scfDG.hamDG.density{elemIdx} .* ...
                (sumPseudoCharge / sumDensity);
        end
    end
end     

% ----------------------- end of density -------------------------------


% -------------- wavefunctions in extended element ---------------------

if esdfParam.userOption.general.isRestartWfn
    restartFilePrefix = scfDG.dataFileIO.restartWfn;
    for elemIdx = 1 : numElemTotal
        fileName = restartFilePrefix + "_" + num2str(elemIdx) + ".mat";
        restartData = load(fileName);
        wfnExtElem = restartData.wfnExtElem;

        if size(wfnExtElem) ~= size(scfDG.vecEigSol{elemIdx}.psi.wavefun)
            error('The size of restarting wavefun does not match with current setup');
        end
        scfDG.vecEigSol{elemIdx}.psi.wavefun = wfnExtElem;    
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

dmExtElem = scfDG.vecEigSol{1}.fft.domain;
dmElem = Domain();
for d = 1 : dimDef()
    dmElem.length(d) = scfDG.domain.length(d) / numElem(d);
    dmElem.numGrid(d) = scfDG.domain.numGrid(d) / numElem(d);
    dmElem.numGridFine(d) = scfDG.domain.numGridFine(d) / numElem(d);
    % posStart relative to the extended element
    dmExtElem.posStart(d) = 0;
    if numElem(d) > 1
        dmElem.posStart(d) = dmElem.length(d) * esdfParam.DG.bufferSize;
    else
        dmElem.posStart(d) = 0;
    end
end

numLGL             = hamDG.grid.numLGLGridElem;
lengthUniform      = dmExtElem.length;
    
LGLGrid             = LGLMesh(dmElem, numLGL);
UniformGrid         = UniformMesh(dmExtElem);
UniformGridFine     = UniformMeshFine(dmExtElem);

scfDG.PeriodicUniformToLGLMat = GenerateFourierInterpMat(...
    LGLGrid, UniformGrid, lengthUniform);
scfDG.PeriodicUniformFineToLGLMat = GenerateFourierInterpMat(...
    LGLGrid, UniformGridFine, lengthUniform);

% ----------------- end of grid transfer matrix -------------------------


% --------------------- periodize the potential ------------------------

% Whether to periodize the potential in the extended element.
if esdfParam.userOption.DG.isPeriodizePotential
    DIM = dimDef();
    scfDG.periodicPotential.vBubble = cell(DIM, 1);
    distancePeriodize = scfDG.periodicPotential.distancePeriodize;
    
    dmExtElem = scfDG.vecEigSol{1}.fft.domain;
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