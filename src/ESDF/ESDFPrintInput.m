function ESDFPrintInput(esdfParam)
% ESDFPRINTINPUT prints ESDF parameters in ESDFInputParam object esdfParam 
%    to output file with ID outFid, if outFid == 1, print to the screen.
%
%    See also ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% If the product of the number of elements is 1, recognize this as a PWDFT
% calculation.

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
InfoPrint(outFid, "Number of Extra Electron             = ",  esdfParam.basic.extraElectron);

    
% DGDFT or PWDFT
if esdfParam.isDGDFT
    PrintBlock(outFid, "DGDFT information");
    
    InfoPrint(outFid, "Penalty Alpha                        = ",  esdfParam.DG.penaltyAlpha);
    InfoPrint(outFid, "Element size                         = ",  esdfParam.DG.numElem); 
    InfoPrint(outFid, "Wfn Elem GridSize                    = ",  esdfParam.DG.numGridWavefunctionElem);
    InfoPrint(outFid, "Rho Elem GridSize                    = ",  esdfParam.DG.numGridDensityElem); 
    InfoPrint(outFid, "LGL Grid size                        = ",  esdfParam.DG.numGridLGL); 
    InfoPrint(outFid, "LGL GridFactor                       = ",  esdfParam.DG.LGLGridFactor);
    InfoPrint(outFid, "Buffer size                          = ",  esdfParam.DG.bufferSize);

    InfoPrint(outFid, "SVD Basis Tol                        = ",  esdfParam.control.SVDBasisTolerance);
    InfoPrint(outFid, "SCF Inner Tol                        = ",  esdfParam.control.scfInnerTolerance);
    InfoPrint(outFid, "SCF Inner MaxIter                    = ",  esdfParam.control.scfInnerMaxIter);
    
    InfoPrint(outFid, "DG Solver                            = ",  esdfParam.basic.DGSolver);    
    
    InfoPrint(outFid, "A Posteriori error each step         = ",  esdfParam.userOption.DG.isCalculateAPosterioriEachSCF);

    if esdfParam.userOption.DG.isPeriodizePotential
        InfoPrint(outFid, "PeriodizePotential                   = ",  esdfParam.userOption.DG.isPeriodizePotential);
        InfoPrint(outFid, "DistancePeriodize                    = ",  esdfParam.DG.distancePeriodize);
    end

    InfoPrint(outFid, 'Number of ALB for element: %d \n', esdfParam.DG.numALBElem);
    
else  % PWDFT
    PrintBlock(outFid, "PWDFT information");
    
    if esdfParam.basic.PWSolver == "PPCG" 
        InfoPrint(outFid, "Subblock size (sbSize) in PPCG       = ",  esdfParam.PW.PPCGsbSize); 
    end
    InfoPrint(outFid, "OutputWavefun                        = ",  esdfParam.userOption.PW.isOutputWavefun);
    InfoPrint(outFid, "OutputDensity                        = ",  esdfParam.userOption.PW.isOutputDensity);
    InfoPrint(outFid, "OutputPotential                      = ",  esdfParam.userOption.PW.isOutputPotential);

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