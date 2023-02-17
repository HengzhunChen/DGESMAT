function HamKS = CalculateForce(HamKS, psi)
% HAMILTONNIANKS/CALCULATEFORCE calculates force for each atom.
%
%    HamKS = CalculateForce(HamKS, psi) computes force for each atom from 
%    wavefun psi and stores the result in HamKS.atomList(:).force
%    
%    See also HamiltonianKS, Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


dim = dimDef();
vol = HamKS.domain.Volume();
numGridTotalFine = HamKS.domain.NumGridTotalFine();

numState = psi.NumStateTotal();
numAtom = length(HamKS.atomList);

force = zeros(numAtom, dim);
F = HamKS.fft;


% *********************************************************************
% Compute the force from local pseudopotential
% *********************************************************************

% Using integration by parts for local pseudopotential.
% No need to evaluate the derivative of the local pseudopotential.
% This could potentially save some coding effort and perhaps better for
% other pseudopotential such as Troullier-Martins

if ~HamKS.userOption.isUseVLocal
    % pseudocharge formulation of the local contribution to the force
    vhartDrv = cell(dim, 1);
        
    totalCharge = HamKS.density - HamKS.pseudoCharge;
    
    % total charge in the Fourier space    
    totalChargeFourier = F * totalCharge;
    
    % Compute the derivative of the Hartree potential via Fourier transform
    gkkFine = F.gkkFine;
    idxnz = gkkFine ~= 0;
    for d = 1 : dim
        ikFine = F.ikFine{d};
        
        y = totalChargeFourier;
        y(~idxnz) = 0;
        y(idxnz) = y(idxnz) .* 4*pi ./ gkkFine(idxnz) .* ikFine(idxnz);        
        x = F' * y;
        
        % vhartDrv saves the derivative of the Hartree potential
        vhartDrv{d} = real(x);
    end

    for i = 1 : numAtom
        pp = HamKS.pseudoList(i);
        sp = pp.pseudoCharge;
        idx = sp.idx;
        val = sp.val;
        
        wgt = vol / numGridTotalFine;
        resX = sum( val .* vhartDrv{1}(idx) ) * wgt;
        resY = sum( val .* vhartDrv{2}(idx) ) * wgt;
        resZ = sum( val .* vhartDrv{3}(idx) ) * wgt;
        
        force(i, 1) = force(i, 1) + resX;
        force(i, 2) = force(i, 2) + resY;
        force(i, 3) = force(i, 3) + resZ;
    end
    
else
    % VLocal formulation of the local contribution to the force
    
    % first contribution from the pseudocharge
    vhartDrv = cell(dim, 1);
    
    totalCharge = HamKS.density - HamKS.pseudoCharge;
        
    % total charge in the Fourier space
    totalChargeFourier = F * totalCharge;
    
    % compute the derivative of the Hartree potential via Fourier transform
    gkkFine = F.gkkFine;
    idxnz = gkkFine ~= 0;
    for d = 1 : dim
        ikFine = F.ikFine{d};
        y = totalChargeFourier;
        y(~idxnz) = 0;
        y(idxnz) = y(idxnz) .* 4*pi ./ gkkFine(idxnz) .* ikFine(idxnz);
        x = F' * y;
        
        % vhartDrv saves the derivative of the Hartree potential
        vhartDrv{d} = real(x);
    end

    for i = 1 : numAtom
        pp = HamKS.pseudoList(i);
        sp = pp.pseudoCharge;
        idx = sp.idx;
        val = sp.val;
        
        wgt = vol / numGridTotalFine;
        resX = sum( val .* vhartDrv{1}(idx) ) * wgt;
        resY = sum( val .* vhartDrv{2}(idx) ) * wgt;
        resZ = sum( val .* vhartDrv{3}(idx) ) * wgt;
        
        force(i, 1) = force(i, 1) + resX;
        force(i, 2) = force(i, 2) + resY;
        force(i, 3) = force(i, 3) + resZ;
    end
    
    % Second, contribution from vLocalSR   
    % The integration by parts formula requires the grad density
    HamKS.gradDensity = CalculateGradDensity(HamKS);
    
    for i = 1 : numAtom
        pp = HamKS.pseudoList(i);
        sp = pp.vLocalSR;
        idx = sp.idx;
        val = sp.val;
        
        wgt = vol / numGridTotalFine;
        resX = -sum( val .* HamKS.gradDensity{d}(idx) ) * wgt;
        resY = -sum( val .* HamKS.gradDensity{d}(idx) ) * wgt;
        resZ = -sum( val .* HamKS.gradDensity{d}(idx) ) * wgt;
        
        force(i, 1) = force(i, 1) + resX;
        force(i, 2) = force(i, 2) + resY;
        force(i, 3) = force(i, 3) + resZ;
    end
    
end


% *********************************************************************
% Compute the force from nonlocal pseudopotential
% *********************************************************************

% FIXME
% Using intergration by parts, and throw the derivative to the
% wavefunctions.
% No need to evaluate the derivative of the non-local pseudopotential.
% This could potentially save some coding effort, and perhaps better for
% other pseudopotential such as Troullier-Martins

psiDrvFine = cell(dim, 1);

% Loop over atoms and pseudopotentials
for g = 1 : numState
    % compute the derivative of the wavefunctions on a fine grid

    psiFine = CoarseToFine(psi(:, g), F);
    psiFine = psiFine.wavefun;
    
    psiFourier = F * psiFine;
    
    % derivative of psi on a fine grid
    for d = 1 : dim
        ikFine = F.ikFine{d};
        y = psiFourier;
        x = F' * (y .* ikFine);
        psiDrvFine{d} = real(x);
    end


    % Evaluate the contribution to the atomic force
    wgt = vol / numGridTotalFine;
    occupationRate = HamKS.occupationRate(g);

    for i = 1 : numAtom
        vnlList = HamKS.pseudoList(i).vnlList;
        if ~isempty(vnlList)
            idx = vnlList.idx;
            val = vnlList.val;
            gamma = vnlList.wgt;
    
            resVAL = val' * psiFine(idx);
            resDX = val' * psiDrvFine{1}(idx);
            resDY = val' * psiDrvFine{2}(idx);
            resDZ = val' * psiDrvFine{3}(idx);
    
            force(i, 1) = force(i, 1) - ...
                4 * occupationRate * wgt * sum(gamma .* resVAL .* resDX);
            force(i, 2) = force(i, 2) - ...
                4 * occupationRate * wgt * sum(gamma .* resVAL .* resDY);
            force(i, 3) = force(i, 3) - ...
                4 * occupationRate * wgt * sum(gamma .* resVAL .* resDZ);
        end
    end  % end for i
    
end  % end for g


% *********************************************************************
% Compute the total force and give the value to atomList
% *********************************************************************

for i = 1 : numAtom
    HamKS.atomList(i).force = force(i, :);
end

% add extra contribution to the force
if HamKS.VDWType == "DFT-D2"
    % update force
    for i = 1 : numAtom
        HamKS.atomList(i).force = HamKS.atomList(i).force + HamKS.forceVdw(i, :);
    end
end

% add extra contribution from short range interaction
if HamKS.userOption.isUseVLocal
    for i = 1 : numAtom
        HamKS.atomList(i).force = HamKS.atomList(i).force + HamKS.forceIonSR(i, :);
    end
end

% add the contribution from external force
for i = 1 : numAtom
    HamKS.atomList(i).force = HamKS.atomList(i).force + HamKS.forceExt(i, :);
end   

end