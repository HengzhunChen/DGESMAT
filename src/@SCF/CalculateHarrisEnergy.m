function EfreeHarris = CalculateHarrisEnergy(scf)
% SCF/CALCULATEHARRISENERGY computes Harris energy.
%
%    NOTE:
%    The difference between the Kohn-Sham energy and the Harris energy
%    is that the nonlinear correction term in the Harris energy
%    functional must be computed via the input electron density, rather
%    than the output electron density or the mixed electron density.
%   
%    Reference:
%   
%    [Soler et al. "The SIESTA method for ab initio order-N
%    materials", J. Phys. Condens. Matter. 14, 2745 (2002) pp 18]
%
%    See also SCF.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


hamKS = scf.eigSol.hamKS;

eigVal = hamKS.eigVal;
occupationRate = hamKS.occupationRate;

% Kinetic energy %
numSpin = hamKS.numSpin;
Ekin = numSpin * sum(eigVal .* occupationRate);

% Self energy %
Eself = hamKS.Eself;

% Ionic repulsion related energy %
EIonSR = hamKS.EIonSR;

% Van der Waals energy
EVdw = hamKS.EVdw;

% External energy
Eext = hamKS.Eext;

% Nonlinear correction part. This part used the Hartree energy and XC
% correlation energy from the old electron density.
density = hamKS.density;
vxc = hamKS.vxc;
pseudoCharge = hamKS.pseudoCharge;
vhart = hamKS.vhart;

EVxc = sum( vxc .* density );
Ehart = 0.5 * sum( vhart .* (density + pseudoCharge) );

ntotFine = hamKS.domain.NumGridTotalFine();
vol = hamKS.fft.domain.Volume();

Ehart = Ehart * vol / ntotFine;
EVxc = EVxc * vol / ntotFine;

Exc = scf.Exc;

% Correction energy
Ecor = (Exc - EVxc) - Ehart - Eself + EIonSR + EVdw + Eext;

% Helmholtz free energy
if hamKS.numOccupiedState == hamKS.NumStateTotal()
    % zero temperature
    Efree = Ekin + Ecor;
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
    Efree = Efree + Ecor + ...
        fermi * hamKS.numOccupiedState * numSpin;
end

EfreeHarris = Efree;

end