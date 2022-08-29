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

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% TODO: the code for smearing have NOT be tested

hamDG = scfDG.hamDG;
eigVal = hamDG.eigVal;
occupationRate = hamDG.occupationRate;

numElem = scfDG.numElem;

% NOTE: To avoid confusion, all energies in this routine are temporary
% variables other than scfDG.EfreeHarris.
%
% The related energies will be computed again in the routine 
% CalculateKSEnergy()

% kinetic energy from the new density matrix
numSpin = hamDG.numSpin;

if scfDG.CheFSIDG.isUseCompSubspace    
    %
    % TODO
    %
else
    Ekin = numSpin * sum(eigVal .* occupationRate);    
end


% self energy part
Eself = 0;
atomList = hamDG.atomList;
for i = 1 : length(atomList)
    type = atomList(i).type;
    Eself = Eself + scfDG.ptable.SelfIonInteraction(type);
end


% Nonlieanr correction part. This part uses the Hartree energy and XC
% correlation energy from the old electron density.

Ehart = 0;
EVxc = 0;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            density      = hamDG.density{i, j, k};
            vxc          = hamDG.vxc{i, j, k};
            pseudoCharge = hamDG.pseudoCharge{i, j, k};
            vhart        = hamDG.vhart{i, j, k};
            
            EVxc = EVxc + sum(vxc .* density);
            Ehart = Ehart + 0.5 * sum( vhart .* (density + pseudoCharge) );
        end
    end
end

Ehart = Ehart * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
EVxc  = EVxc  * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();

% use the previous exchange-correlation energy
Exc = scfDG.Exc;

% Correction energy
Ecor = (Exc - EVxc) - Ehart - Eself;

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
    
    if scfDG.CheFSIDG.isUseCompSubspace
        % Complementary subspace technique in use
        %
        % TODO
        %
    else
        % Complementary subspace technique not in use : full spectrum
        % available
        if SmearingScheme == "FD"
            idx = eigVal >= fermi;
            EfreeHarris = EfreeHarris - numSpin / Tbeta * ...
                                sum( log(1 + exp(-Tbeta*(eigVal(idx) - fermi))) );
            EfreeHarris = EfreeHarris + ...
                                numSpin * sum( (eigVal(~idx) - fermi) ) -...
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
        
    end  % end of full spectrum available calculation
end  % end of finite temperature calculation
            

end