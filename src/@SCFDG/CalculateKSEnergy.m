function scfDG = CalculateKSEnergy(scfDG)
% SCFDG/CALCULATEKSENERGY calculate Kohn-Sham DFT energies.
%
%    scfDG = CalculateKSEnergy(scfDG) computes most DFT related energies
%    and save them back to scfDG.
%
%    See also SCFDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% TODO: the code for Efree have NOT be tested

hamDG = scfDG.hamDG;
eigVal = hamDG.eigVal;
occupationRate = hamDG.occupationRate;


% Band energy
numSpin = hamDG.numSpin;
Ekin = numSpin * sum(eigVal .* occupationRate);
scfDG.Ekin = Ekin;

% Hartree and XC part
numElem = scfDG.numElem;
numElemTotal = prod(numElem);

Ehart = 0;
EVxc = 0;

for elemIdx = 1 : numElemTotal
    density      = scfDG.hamDG.density{elemIdx};
    vxc          = scfDG.hamDG.vxc{elemIdx};
    pseudoCharge = scfDG.hamDG.pseudoCharge{elemIdx};
    vhart        = scfDG.hamDG.vhart{elemIdx};
    
    EVxc = EVxc + sum(vxc .* density);
    % NOTE the sign flip
    Ehart = Ehart + 0.5 * sum( vhart .* (density + pseudoCharge) );
end

Ehart = Ehart * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
EVxc  = EVxc  * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();

scfDG.Ehart = Ehart;
scfDG.EVxc = EVxc;

% Self energy and ionic repulsion energy
scfDG.Eself = hamDG.Eself;
scfDG.EIonSR = hamDG.EIonSR;

% Van der waals energy
scfDG.Evdw = hamDG.EVdw;

% External energy
scfDG.Eext = hamDG.Eext;

% Correction energy
Ecor = (scfDG.Exc - EVxc) - Ehart - scfDG.Eself + scfDG.EIonSR + scfDG.Evdw + scfDG.Eext;
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