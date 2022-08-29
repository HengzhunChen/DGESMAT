function ESDFPrintInput()
% ESDFPRINTINPUT prints ESDF parameters in global variable esdfParam to
%    output file with ID outFid, if outFid == 1, print to the screen.
%
%    See also ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% If the product of the number of elements is 1, recognize this as a PWDFT
% calculation.

global esdfParam;

outFid = 0;

PrintBlock(outFid, "Common information");

InfoPrint(outFid, "Super cell                           = ",  esdfParam.basic.domain.length);
InfoPrint(outFid, "Grid Wavefunction                    = ",  esdfParam.basic.domain.numGrid); 
InfoPrint(outFid, "Grid Density                         = ",  esdfParam.basic.domain.numGridFine);
InfoPrint(outFid, "Mixing dimension                     = ",  esdfParam.basic.mixMaxDim);
InfoPrint(outFid, "Mixing variable                      = ",  esdfParam.basic.mixVariable);
InfoPrint(outFid, "Mixing type                          = ",  esdfParam.basic.mixType);
InfoPrint(outFid, "Mixing Steplength                    = ",  esdfParam.basic.mixStepLength);

InfoPrint(outFid, "Temperature                          = ",  au2KDef() / esdfParam.basic.Tbeta, "[K]");
InfoPrint(outFid, "Extra states                         = ",  esdfParam.basic.numExtraState);
InfoPrint(outFid, "Smearing scheme                      = ",  esdfParam.basic.smearingScheme);
InfoPrint(outFid, "PW Solver                            = ",  esdfParam.basic.PWSolver);
InfoPrint(outFid, "XC Type                              = ",  esdfParam.basic.XCType);
InfoPrint(outFid, "Pseudo Type                          = ",  esdfParam.basic.pseudoType);
for i = 1 : length(esdfParam.basic.upfFile)
    InfoPrint(outFid, "UPF File                             = ",  esdfParam.basic.upfFile{i});
end

InfoPrint(outFid, "SCF Outer Tol                        = ",  esdfParam.control.scfOuterTolerance);
InfoPrint(outFid, "SCF Outer MaxIter                    = ",  esdfParam.control.scfOuterMaxIter);
InfoPrint(outFid, "SCF Free Energy Per Atom Tol         = ",  esdfParam.control.scfOuterEnergyTolerance);
InfoPrint(outFid, "Eig Min Tolerence                    = ",  esdfParam.control.eigMinTolerance);
InfoPrint(outFid, "Eig Tolerence                        = ",  esdfParam.control.eigTolerance);
InfoPrint(outFid, "Eig MaxIter                          = ",  esdfParam.control.eigMaxIter);
InfoPrint(outFid, "Eig Tolerance Dyn                    = ",  esdfParam.userOption.general.isPWeigTolDynamic);
InfoPrint(outFid, "Num unused state                     = ",  esdfParam.basic.numUnusedState);
InfoPrint(outFid, "EcutWavefunction                     = ",  esdfParam.basic.ecutWavefunction);
InfoPrint(outFid, "Density GridFactor                   = ",  esdfParam.basic.densityGridFactor);

InfoPrint(outFid, "Use Atom Density                     = ",  esdfParam.userOption.general.isUseAtomDensity);
InfoPrint(outFid, "Use VLocal                           = ",  esdfParam.userOption.general.isUseVLocal);
InfoPrint(outFid, "RestartDensity                       = ",  esdfParam.userOption.general.isRestartDensity);
InfoPrint(outFid, "RestartWfn                           = ",  esdfParam.userOption.general.isRestartWfn);
InfoPrint(outFid, "OutputDensity                        = ",  esdfParam.userOption.general.isOutputDensity);
InfoPrint(outFid, "Number of Extra Electron             = ",  esdfParam.basic.extraElectron);


% Ionic motion
if ~isempty(esdfParam.ionDyn.ionMove)
    InfoPrint(outFid, "");
    InfoPrint(outFid, "Ion move mode                        = ",  esdfParam.ionDyn.ionMove);
    InfoPrint(outFid, "Max steps for ion                    = ",  esdfParam.ionDyn.ionMaxIter);
    InfoPrint(outFid, "MD time step                         = ",  esdfParam.ionDyn.MDTimeStep);
    InfoPrint(outFid, "Ion Temperature                      = ",  esdfParam.ionDyn.ionTemperature, "[K]");
    InfoPrint(outFid, "Thermostat mass                      = ",  esdfParam.ionDyn.qMass);
    InfoPrint(outFid, "Langevin damping                     = ",  esdfParam.ionDyn.langevinDamping);
    InfoPrint(outFid, "RestartPosition                      = ",  esdfParam.userOption.ionDyn.isRestartPosition);
    InfoPrint(outFid, "RestartVelocity                      = ",  esdfParam.userOption.ionDyn.isRestartVelocity);
    InfoPrint(outFid, "OutputPosition                       = ",  esdfParam.userOption.ionDyn.isOutputPosition);
    InfoPrint(outFid, "OutputVelocity                       = ",  esdfParam.userOption.ionDyn.isOutputVelocity);
    InfoPrint(outFid, "Output XYZ format                    = ",  esdfParam.userOption.ionDyn.isOutputXYZ);
    InfoPrint(outFid, "Force tol for geoopt                 = ",  esdfParam.ionDyn.geoOptMaxForce);
    InfoPrint(outFid, "MD extrapolation type                = ",  esdfParam.ionDyn.MDExtrapolationType);
    InfoPrint(outFid, "MD extrapolation variable            = ",  esdfParam.ionDyn.MDExtrapolationVariable);
    InfoPrint(outFid, "MD SCF Phi MaxIter                   = ",  esdfParam.hybrid.MDscfPhiMaxIter);
    InfoPrint(outFid, "MD SCF Outer MaxIter                 = ",  esdfParam.hybrid.MDscfOuterMaxIter);
    InfoPrint(outFid, "MD SCF Energy Criteria Engage Iter   = ",  esdfParam.ionDyn.MDscfEnergyCriteriaEngageIonIter);
    InfoPrint(outFid, "MD SCF Etot diff                     = ",  esdfParam.ionDyn.MDscfEtotDiffTol);
    InfoPrint(outFid, "MD SCF Eband diff                    = ",  esdfParam.ionDyn.MDscfEbandDiffTol);
    InfoPrint(outFid, "");
end
    
% DGDFT or PWDFT
if esdfParam.isDGDFT
    PrintBlock(outFid, "DGDFT information");
    
    InfoPrint(outFid, "Penalty Alpha                        = ",  esdfParam.DG.penaltyAlpha);
    InfoPrint(outFid, "Element size                         = ",  esdfParam.DG.numElem); 
    InfoPrint(outFid, "Wfn Elem GridSize                    = ",  esdfParam.DG.numGridWavefunctionElem);
    InfoPrint(outFid, "Rho Elem GridSize                    = ",  esdfParam.DG.numGridDensityElem); 
    InfoPrint(outFid, "LGL Grid size                        = ",  esdfParam.DG.numGridLGL); 
    InfoPrint(outFid, "LGL GridFactor                       = ",  esdfParam.DG.LGLGridFactor);

    InfoPrint(outFid, "SVD Basis Tol                        = ",  esdfParam.control.SVDBasisTolerance);
    InfoPrint(outFid, "SCF Inner Tol                        = ",  esdfParam.control.scfInnerTolerance);
    InfoPrint(outFid, "SCF Inner MaxIter                    = ",  esdfParam.control.scfInnerMaxIter);
    
    InfoPrint(outFid, "DG Solver                            = ",  esdfParam.basic.DGSolver);    
    
    InfoPrint(outFid, "OutputWfnExtElem                     = ",  esdfParam.userOption.DG.isOutputWfnExtElem);
    InfoPrint(outFid, "Force each step (deprecated)         = ",  esdfParam.userOption.DG.isCalculateForceEachSCF); 
    InfoPrint(outFid, "A Posteriori error each step         = ",  esdfParam.userOption.DG.isCalculateAPosterioriEachSCF);

    if esdfParam.userOption.DG.isPeriodizePotential
        InfoPrint(outFid, "PeriodizePotential                 = ",  esdfParam.userOption.DG.isPeriodizePotential);
        InfoPrint(outFid, "DistancePeriodize                  = ",  esdfParam.DG.distancePeriodize);
    end
    % FIXME Potentially obsolete potential barriers
    if esdfParam.userOption.DG.isPotentialBarrier
        InfoPrint(outFid, "Potential Barrier = ",  esdfParam.userOption.DG.isPotentialBarrier);
        InfoPrint(outFid, "Barrier W         = ",  esdfParam.DG.potentialBarrierW);
        InfoPrint(outFid, "Barrier S         = ",  esdfParam.DG.potentialBarrierS);
        InfoPrint(outFid, "Barrier R         = ",  esdfParam.DG.potentialBarrierR);
    end

    InfoPrint(outFid, 'Number of ALB for element: %d \n', esdfParam.DG.numALBElem);
    
else  % PWDFT
    PrintBlock(outFid, "PWDFT information");

    % FIXME For DG as well later
    InfoPrint(outFid, "SCF Phi MaxIter                      = ",  esdfParam.hybrid.scfPhiMaxIter);
    InfoPrint(outFid, "SCF Phi Tol                          = ",  esdfParam.hybrid.scfPhiTolerance);
    InfoPrint(outFid, "Hybrid ACE                           = ",  esdfParam.hybrid.isHybridACE);
    InfoPrint(outFid, "Hybrid DF                            = ",  esdfParam.hybrid.isHybridDF);
    InfoPrint(outFid, "Hybrid DF Type                       = ",  esdfParam.hybrid.DFType);
    InfoPrint(outFid, "Hybrid DF Kmeans WF Type             = ",  esdfParam.hybrid.DFKmeansWFType);
    InfoPrint(outFid, "Hybrid DF Kmeans WF Alpha            = ",  esdfParam.hybrid.DFKmeansWFAlpha);
    InfoPrint(outFid, "Hybrid Active Init                   = ",  esdfParam.hybrid.isHybridActiveInit);
    InfoPrint(outFid, "Hybrid Mixing Type                   = ",  esdfParam.hybrid.MixType);
    InfoPrint(outFid, "Hybrid DF Num Mu                     = ",  esdfParam.hybrid.DFNumMu);
    InfoPrint(outFid, "Hybrid DF Num GaussianRandom         = ",  esdfParam.hybrid.DFNumGaussianRandom);
    InfoPrint(outFid, "Hybrid DF Tolerance                  = ",  esdfParam.hybrid.DFTolerance);
    InfoPrint(outFid, "EXX div type                         = ",  esdfParam.hybrid.exxDivergenceType);
    
    if esdfParam.basic.PWSolver == "PPCG" 
        InfoPrint(outFid, "Subblock size (sbSize) in PPCG       = ",  esdfParam.PW.PPCGsbSize); 
    end
    InfoPrint(outFid, "OutputWfn                            = ",  esdfParam.userOption.PW.isOutputWfn);

end

% Information of atoms
atomList = esdfParam.basic.atomList;
InfoPrint(outFid, ""); 
InfoPrint(outFid, "NumAtom = ", length(atomList)); 
InfoPrint(outFid, "Atom Type and Coordinates (possibly after shifting into unit cell)");
InfoPrint(outFid, ""); 
for i = 1 : length(atomList)
    InfoPrint(outFid, "Type = ", atomList(i).type, "Position = ", atomList(i).pos);
end
InfoPrint(outFid, ""); 

end