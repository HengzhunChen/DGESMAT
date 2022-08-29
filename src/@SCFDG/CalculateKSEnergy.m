function scfDG = CalculateKSEnergy(scfDG)
% SCFDG/CALCULATEKSENERGY calculate Kohn-Sham DFT energies.
%
%    scfDG = CalculateKSEnergy(scfDG) computes most DFT related energies
%    and save them back to scfDG.
%
%    See also SCFDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% TODO: the code for Efree have NOT be tested

hamDG = scfDG.hamDG;
eigVal = hamDG.eigVal;
occupationRate = hamDG.occupationRate;


% Band energy
numSpin = hamDG.numSpin;

if scfDG.CheFSIDG.isUseCompSubspace
    %
    % TODO
    %
else
    Ekin = numSpin * sum(eigVal .* occupationRate);
end
scfDG.Ekin = Ekin;

% Self energy part
Eself = 0;
atomList = hamDG.atomList;
for i = 1 : length(atomList)
    type = atomList(i).type;
    Eself = Eself + scfDG.ptable.SelfIonInteraction(type);
end
scfDG.Eself = Eself;

% Hartree and XC part
numElem = scfDG.numElem;

Ehart = 0;
EVxc = 0;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            density      = scfDG.hamDG.density{i, j, k};
            vxc          = scfDG.hamDG.vxc{i, j, k};
            pseudoCharge = scfDG.hamDG.pseudoCharge{i, j, k};
            vhart        = scfDG.hamDG.vhart{i, j, k};
            
            EVxc = EVxc + sum(vxc .* density);
            % NOTE the sign flip
            Ehart = Ehart + 0.5 * sum( vhart .* (density + pseudoCharge) );
        end
    end
end

Ehart = Ehart * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
EVxc  = EVxc  * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();

scfDG.Ehart = Ehart;
scfDG.EVxc = EVxc;

% Correction energy
Ecor = (scfDG.Exc - EVxc) - Ehart - Eself;
scfDG.Ecor = Ecor;

% Total energy
scfDG.Etot = scfDG.Ekin + scfDG.Ecor;

% Helmholtz free energy
if hamDG.numOccupiedState == hamDG.NumStateTotal()
    % zero temperature
    scfDG.Efree = scfDG.Etot;
else
    % finite temperature
    scfDG.Efree = 0;
    fermi = scfDG.fermi;
    Tbeta = scfDG.Tbeta;
    
    if scfDG.CheFSIDG.isUseCompSubspace
        %
        % TODO
        %
    else
        idx = eigVal >= fermi;
        scfDG.Efree = scfDG.Efree - numSpin / Tbeta * ...
                      sum( log(1 + exp(-Tbeta*(eigVal(idx) - fermi))) );
        scfDG.Efree = scfDG.Efree + ...
                      numSpin * sum( (eigVal(~idx) - fermi) ) - ...
                      numSpin / Tbeta * sum( log(1 + exp(Tbeta*(eigVal(~idx) - fermi))) );

        scfDG.Efree = scfDG.Efree + Ecor + ...
                      fermi * hamDG.numOccupiedState * numSpin;
    end
end


end