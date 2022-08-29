function PrintState(scf, iter)
% SCF/PRINTSTATE print state information for iter-th SCF iteration. 
%
%    See also SCF.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


HOMO = scf.eigSol.eigVal( scf.eigSol.hamKS.numOccupiedState );
numExtraState = scf.eigSol.hamKS.numExtraState;
if numExtraState > 0
    % lowest unoccupied molecular orbital (LUMO)
    LUMO = scf.eigSol.eigVal( scf.eigSol.hamKS.numOccupiedState+1 );
end

for i = 1 : length(scf.eigSol.eigVal)
    InfoPrint(0, 'band#  = ', i, ...
                 'eigval  = ', scf.eigSol.eigVal(i), ...
                 'resval  = ', scf.eigSol.resVal(i), ...
                 'occrate  = ', scf.eigSol.hamKS.occupationRate(i) );
end

InfoPrint(0, '\n');
InfoPrint(0, 'NOTE:  Ecor  = Exc - EVxc - Ehart -Eself + EIonSR + EVdw + Eext \n');
InfoPrint(0, '       Etot  = Ekin + Ecor \n');
InfoPrint(0, '       Efree = Etot + Entropy \n \n');

InfoPrint([0, 1], 'SCF Iteration # %d \n', iter);

InfoPrint([0, 1], "Etot              = ",  scf.Etot, "[au]");
InfoPrint(0, "Efree             = ",  scf.Efree, "[au]");
InfoPrint(0, "EfreeHarris       = ",  scf.EfreeHarris, "[au]");
InfoPrint(0, "Ekin              = ",  scf.Ekin, "[au]");
InfoPrint(0, "Ehart             = ",  scf.Ehart, "[au]");
InfoPrint(0, "EVxc              = ",  scf.EVxc, "[au]");
InfoPrint(0, "Exc               = ",  scf.Exc, "[au]"); 
InfoPrint(0, "EVdw              = ",  scf.EVdw, "[au]"); 
InfoPrint(0, "Eself             = ",  scf.Eself, "[au]");
InfoPrint(0, "EIonSR            = ",  scf.EIonSR, "[au]");
InfoPrint(0, "Eext              = ",  scf.Eext, "[au]");
InfoPrint(0, "Ecor              = ",  scf.Ecor, "[au]");
InfoPrint(0, "Fermi             = ",  scf.fermi, "[au]");
InfoPrint(0, "Total charge      = ",  scf.totalCharge, "[au]");
InfoPrint(0, "HOMO              = ",  HOMO*au2evDef(), "[eV]");

if numExtraState > 0 
    InfoPrint(0, "LUMO              = ", LUMO*au2evDef(), "[eV]");
end

InfoPrint(0, '\n');

end