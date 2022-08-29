function scf = CalculateEnergy(scf)
% SCF/CALCULATEENERGY computes most energies in PWDFT and save them in SCF 
%    object scf.
%
%    See also SCF.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam

hamKS = scf.eigSol.hamKS;

eigVal = hamKS.eigVal;
numSpin = hamKS.numSpin;

% Kinetic energy %
occupationRate = hamKS.occupationRate;
scf.Ekin = numSpin * sum(eigVal .* occupationRate);

% Hartree and xc part %
ntot = hamKS.domain.NumGridTotalFine();
vol = hamKS.domain.Volume();

density = hamKS.density;
vxc = hamKS.vxc;
pseudoCharge = hamKS.pseudoCharge;
vhart = hamKS.vhart;

EVxc = sum(vxc .* density);
Ehart = 0.5 * sum( vhart .* (density + pseudoCharge) );

scf.Ehart = Ehart * vol / ntot;
scf.EVxc = EVxc * vol / ntot;

% Ionic repulsion related energy %
scf.Eself = hamKS.Eself;
scf.Ecor = (scf.Exc - scf.EVxc) - scf.Ehart - scf.Eself;
if esdfParam.userOption.general.isUseVLocal
    scf.EIonSR = hamKS.EIonSR;
    scf.Ecor = scf.Ecor +  scf.EIonSR;
end

% Van der waals energy %
scf.EVdw = hamKS.EVdw;
scf.Ecor = scf.Ecor + scf.EVdw;

% External energy %
scf.Eext = hamKS.Eext;
scf.Ecor = scf.Ecor + scf.Eext;

% Total energy %
scf.Etot = scf.Ekin + scf.Ecor;


% Helmholtz free energy %
if hamKS.numOccupiedState == hamKS.NumStateTotal()
    % zero temperature
    scf.Efree = scf.Etot;
else
    % finite temperature
    Efree = 0;
    fermi = scf.fermi;
    Tbeta = scf.Tbeta;
    idx = eigVal >= fermi;
    Efree = Efree - numSpin / Tbeta * ...
        sum( log(1 + exp(-Tbeta*(eigVal(idx) - fermi))) );
    Efree = Efree + numSpin * sum( (eigVal(~idx) - fermi) ) -...
        numSpin / Tbeta * sum( log(1 + exp(Tbeta*(eigVal(~idx) - fermi))) );
    scf.Efree = Efree + scf.Ecor + ...
        fermi * hamKS.numOccupiedState * numSpin;
end

end