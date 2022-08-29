function EfreeSecondOrder = CalculateSecondOrderEnergy(scfDG)
% SCFDG/CALCULATESECONDORDERENERGY calculate second order energy
%
%    EfreeSecondOrder = CalculateSecondOrderEnergy(scfDG) calculates the 
%    second order accurate energy that is applicable to both density and 
%    potential mixing.
%
%    See also SCFDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

% TODO: the code for EfreeSecondOrder have NOT be tested

hamDG = scfDG.hamDG;
eigVal = hamDG.eigVal;
occupationRate = hamDG.occupationRate;

% NOTE: To avoid confusion, all energies in this routine are temporary
% variables other than scfDG.EfreeSecondOrder
%
% This is similar to the situation in CalculateHarrisEnergy()

% Kinetic energy from the new density matrix
numSpin = hamDG.numSpin;

if scfDG.CheFSIDG.isUseCompSubspace
    %
    % TODO
    %
else
    Ekin = numSpin * sum(eigVal .* occupationRate);
end

% Self energy part
Eself = 0;
atomList = hamDG.atomList;
for i = 1 : length(atomList)
    type = atomList(i).type;
    Eself = Eself + scfDG.ptable.SelfIonInteraction(type);
end

% Nonlinear correction part. This part uses the Hartree energy and XC
% correlation energy from the OUTPUT electron density, but the total
% potential is the INPUT one used in the diagonalization process.
% The density is also the OUTPUT density.
%
% NOTE the sign flip in Ehart, which is different from those in KS energy
% functional and Harris energy functional.

numElem = scfDG.numElem;

Ehart = 0;
EVtot = 0;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            density      = scfDG.hamDG.density{i, j, k};
            vext         = scfDG.hamDG.vext{i, j, k};
            vtot         = scfDG.hamDG.vtot{i, j, k};
            pseudoCharge = scfDG.hamDG.pseudoCharge{i, j, k};
            vhart        = scfDG.hamDG.vhart{i, j, k};
            
            EVtot = EVtot + sum( (vtot - vext) .* density );
            % NOTE the sign flip
            Ehart = Ehart + 0.5 * sum( vhart .* (density - pseudoCharge) );
        end
    end
end

Ehart = Ehart * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
EVtot = EVtot * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();

% Use the exchange-correlation energy with respect to the new electron
% density
Exc = scfDG.Exc;

% Correction energy
% NOTE the correction energy in the second order method means differently
% from that in Harris energy functional or the KS energy functional.
Ecor = (Exc + Ehart - Eself) - EVtot;

% Second order accurate free energy functional
if hamDG.numOccupiedState == hamDG.NumStateTotal()
    % zero temperature
    EfreeSecondOrder = Ekin + Ecor;
else
    % finite temperature
    EfreeSecondOrder = 0;
    fermi = scfDG.fermi;
    Tbeta = scfDG.Tbeta;
    Tsigma = scfDG.Tsigma;

    SmearingScheme = scfDG.smearing.SmearingScheme;
    MPSmearingOrder = scfDG.smearing.MPsmearingOrder;
    
    if scfDG.CheFSIDG.isUseCompSubspace
        %
        % TODO
        %
    else
        % complementary subspace technique not in use: full spectrum
        % available
        if SmearingScheme == "FD"
            idx = eigVal >= fermi;
            EfreeSecondOrder = EfreeSecondOrder - numSpin / Tbeta * ...
                                     sum( log(1 + exp(-Tbeta*(eigVal(idx) - fermi))) );
            EfreeSecondOrder = EfreeSecondOrder + ...
                                numSpin * sum( (eigVal(~idx) - fermi) ) - ...
                                numSpin / Tbeta * sum( log(1 + exp(Tbeta*(eigVal(~idx) - fermi))) );
            
            EfreeSecondOrder = EfreeSecondOrder + Ecor + ...
                                fermi * hamDG.numOccupiedState * numSpin;
        else
            % GB or MP schemes in use            
            occupTol = 1e-12;
            
            occup = occupationRate;
            idx = occup > occupTol && (1.0 - occup) > occupTol;
            
            x = (eigVal(idx) - fermi) / Tsigma;
            occupEnergyPart = sum( MPentropy(x, MPSmearingOrder) );
           
            EfreeSecondOrder = Ekin + Ecor + (numSpin * Tsigma) * occupEnergyPart;
        end
    end  % end of full spectrum available calculation
end  % end of finite temperature calculation

    
end