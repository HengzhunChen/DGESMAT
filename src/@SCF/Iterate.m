function scf = Iterate(scf)
% SCF/ITERATE main iteration of PWDFT.
%
%    See also SCF.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


controlVar = scf.controlVar;

% **********************************************************************
%                            Preparation                               
% **********************************************************************

hamKS = scf.eigSol.hamKS;

% exchange-correlation
if controlVar.isCalculateGradRho
    hamKS.gradDensity = CalculateGradDensity(hamKS);
end
[scf.Exc, hamKS.epsxc, hamKS.vxc] = CalculateXC(hamKS);
% Hartree energy 
hamKS.vhart = CalculateHartree(hamKS);
% total potential
hamKS.vtot = CalculateVtot(hamKS);

scf.eigSol.hamKS = hamKS;


% **********************************************************************
%                            Main Iteration                              
% **********************************************************************
    
% -------------- non-hybrid functional calculation -----------------

PrintBlock(0, 'Starting regular SCF iteration.');
isSCFConverged = false;

for iter = 1 : controlVar.scfMaxIter
    if isSCFConverged
        break;
    end
    
    %***********************************************************
    % Performing each iteration
    %***********************************************************
    
    PrintBlock(0, 'SCF iteration # ', iter);
    
    timeIterStart = tic;
    
    % --------------- solve eigenvalue problem -------------------
    
    % update density, gradDensity, potential (stored in scf.vtotNew)
    scf = InnerSolve(scf, iter);
    
    vtotOld = scf.eigSol.hamKS.vtot;
    normVtotOld = norm(vtotOld, 2);
    normVtotDiff = norm(vtotOld - scf.vtotNew, 2);
    scf.scfNorm = normVtotDiff / normVtotOld;
            
    scf = scf.CalculateEnergy();
    
    scf.PrintState(iter);
    
    numAtom = length(scf.eigSol.hamKS.atomList);
    scf.efreeDifPerAtom = abs(scf.Efree - scf.EfreeHarris) / numAtom;
    
    InfoPrint([0, 1], "norm(out-in)/norm(in) = ", scf.scfNorm);
    InfoPrint([0, 1], "Efree diff per atom   = ", scf.efreeDifPerAtom); 
    
    if scf.scfNorm < controlVar.scfTolerance
        % converged %
        InfoPrint([0, 1], 'SCF is converged in %d steps !\n', iter);
        isSCFConverged = true;
    end
    
    % ----------------- Potential Mixing ----------------------------
    
    if scf.mixing.mixType == "anderson" || scf.mixing.mixType == "kerker+anderson"
        [scf.eigSol.hamKS.vtot, scf.mixing.dfMat, scf.mixing.dvMat] = ...
            scf.AndersonMix(iter, ...
                            vtotOld, ...
                            scf.vtotNew, ...
                            scf.mixing.dfMat, ...
                            scf.mixing.dvMat);
    else
        error('Invalid mixing type');
    end
    
    
    timeIterEnd = toc( timeIterStart );
    InfoPrint(0, 'Total time for this SCF iteration = %8f [s] \n\n', timeIterEnd);
    
    % ------------------- end of iteration -------------------------
end


% ----------------- hybrid functional calculation ---------------------
% TODO


% ---------------------- Calculate Force ------------------------------

scf.eigSol.hamKS = CalculateForce(scf.eigSol.hamKS, scf.eigSol.psi);



% **********************************************************************
%                      Output Information                              
% **********************************************************************

% ----------------------- Print Energy --------------------------------

% highest occupied molecular orbital (HOMO)
HOMO = scf.eigSol.eigVal( scf.eigSol.hamKS.numOccupiedState );
numExtraState = scf.eigSol.hamKS.numExtraState;
if numExtraState > 0
    % lowest unoccupied molecular orbital (LUMO)
    LUMO = scf.eigSol.eigVal( scf.eigSol.hamKS.numOccupiedState+1 );
end

% Print out the energy
PrintBlock([0, 1], 'Energy');

InfoPrint([0, 1], 'NOTE:  Ecor  = Exc - EVxc - Ehart -Eself + EIonSR + EVdw + Eext \n');
InfoPrint([0, 1], '       Etot  = Ekin + Ecor \n');
InfoPrint([0, 1], '       Efree = Etot + Entropy \n \n');

InfoPrint([0, 1], "! Etot              = ",  scf.Etot, "[au]");
InfoPrint([0, 1], "! Efree             = ",  scf.Efree, "[au]");
InfoPrint([0, 1], "! EfreeHarris       = ",  scf.EfreeHarris, "[au]");
InfoPrint([0, 1], "! EVdw              = ",  scf.EVdw, "[au]"); 
InfoPrint([0, 1], "Eext                = ",  scf.Eext, "[au]");
InfoPrint([0, 1], "! Fermi             = ",  scf.fermi, "[au]");
InfoPrint([0, 1], "! HOMO              = ",  HOMO*au2evDef(), "[eV]");

if numExtraState > 0 
    InfoPrint([0, 1], "! LUMO              = ", LUMO*au2evDef(), "[eV]");
end


% --------------------- Print Force ---------------------------------

hamKS = scf.eigSol.hamKS;
PrintBlock([0, 1], 'Atomic Force');
numAtom = length(hamKS.atomList);
forceCM = zeros(1, 3);
for i = 1 : numAtom
    InfoPrint([0, 1], 'atom', i, 'force', hamKS.atomList(i).force);
    forceCM = forceCM + hamKS.atomList(i).force;
end

InfoPrint([0, 1], '\n');
InfoPrint([0, 1], 'force for centroid  : ', forceCM);
InfoPrint([0, 1], 'Max force magnitude : ', MaxForce(hamKS.atomList));
InfoPrint([0, 1], '\n');


% ------------------ Save info and data into files ---------------------

% output structure information as file
if scf.userOption.isOutputAtomStruct
    domain = scf.eigSol.fft.domain;
    atomList = scf.eigSol.hamKS.atomList;
    save(scf.dataFileIO.atomStructPW, 'domain', 'atomList');
end

% output data as files
if scf.userOption.isOutputDensity
    dm = scf.eigSol.fft.domain;
    gridpos = UniformMeshFine(dm);

    % only work for the restricted spin case
    density = scf.eigSol.hamKS.density;
    save(scf.dataFileIO.densityPW, 'gridpos', 'density');
end

if scf.userOption.isOutputWavefun
    dm = scf.eigSol.fft.domain;
    gridpos = UniformMesh(dm);

    wavefun = scf.eigSol.psi.wavefun;
    occupationRate = scf.eigSol.hamKS.occupationRate;
    save(scf.dataFileIO.wavefunPW, 'gridpos', 'wavefun', 'occupationRate');
end

if scf.userOption.isOutputPotential
    dm = scf.eigSol.fft.domain;
    gridpos = UniformMeshFine(dm);
    
    % only work for the restricted spin case
    vtot = scf.eigSol.hamKS.vtot;
    save(scf.dataFileIO.potentialPW, 'gridpos', 'vtot');
end

end