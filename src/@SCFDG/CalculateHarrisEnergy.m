function EfreeHarris = CalculateHarrisEnergy(scfDG)
% SCFDG/CALCULATEHARRISENERGY calculate the Harris (free) energy  
%   
%    EfreeHarris = CalculateHarrisEnergy(scfDG) returns the Harris free
%    energy.
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
%    See also SCFDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% TODO: the code for smearing have NOT be tested

hamDG = scfDG.hamDG;
eigVal = hamDG.eigVal;
occupationRate = hamDG.occupationRate;

numElem = scfDG.numElem;
numElemTotal = prod(numElem);

% NOTE: To avoid confusion, all energies in this routine are temporary
% variables other than scfDG.EfreeHarris.
%
% The related energies will be computed again in the routine 
% CalculateKSEnergy()


% kinetic energy from the new density matrix
numSpin = hamDG.numSpin;
Ekin = numSpin * sum(eigVal .* occupationRate);    

% Nonlieanr correction part. This part uses the Hartree energy and XC
% correlation energy from the old electron density.

Ehart = 0;
EVxc = 0;

for elemIdx = 1 : numElemTotal
    density      = hamDG.density{elemIdx};
    vxc          = hamDG.vxc{elemIdx};
    pseudoCharge = hamDG.pseudoCharge{elemIdx};
    vhart        = hamDG.vhart{elemIdx};
    
    EVxc = EVxc + sum(vxc .* density);
    Ehart = Ehart + 0.5 * sum( vhart .* (density + pseudoCharge) );
end

Ehart = Ehart * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
EVxc  = EVxc  * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();

% use the previous exchange-correlation energy
Exc = scfDG.Exc;

% Self energy and ionic repulsion energy
Eself = hamDG.Eself;
EIonSR = hamDG.EIonSR;

% Van der waals energy
Evdw = hamDG.EVdw;

% External energy
Eext = hamDG.Eext;

% Correction energy
Ecor = (Exc - EVxc) - Ehart - Eself + EIonSR + Evdw + Eext;

% Harris free energy functional
if hamDG.numOccupiedState == hamDG.NumStateTotal()
    % zero temperature
    EfreeHarris = Ekin + Ecor;
else
    % finite temperature
    EfreeHarris = 0;
    
    fermi = scfDG.fermi;
    Tbeta = scfDG.Tbeta;
    Tsigma = scfDG.Tsigma;

    SmearingScheme = scfDG.smearing.SmearingScheme;
    MPsmearingOrder = scfDG.smearing.MPsmearingOrder;
    
    % full spectrum available calculation
    if SmearingScheme == "FD"
        idx = eigVal >= fermi;
        EfreeHarris = EfreeHarris - numSpin / Tbeta * ...
                sum( log(1 + exp(-Tbeta*(eigVal(idx) - fermi))) );
        EfreeHarris = EfreeHarris + ...
                numSpin * sum( (eigVal(~idx) - fermi) ) - ...
                numSpin / Tbeta * sum( log(1 + exp(Tbeta*(eigVal(~idx) - fermi))) );            
        EfreeHarris = EfreeHarris + Ecor + ...
                fermi * hamDG.numOccupiedState * numSpin;
    else
        % GB or MP schemes in use            
        occupTol = 1e-12;
        occup = occupationRate;
        idx = occup > occupTol && (1.0 - occup) > occupTol;
        
        x = (eigVal(idx) - fermi) / Tsigma;
        occupEnergyPart = sum( MPentropy(x, MPsmearingOrder) );
       
        EfreeHarris = Ekin + Ecor + (numSpin * Tsigma) * occupEnergyPart;
    end        
end  % end of finite temperature calculation
            

end