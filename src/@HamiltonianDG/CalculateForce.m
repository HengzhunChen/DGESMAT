function HamDG = CalculateForce(HamDG)
% HAMILTONNIANDG/CALCULATEFORCE calculates force for each atom.
%
%    HamDG = CalculateForce(HamDG) computes force for each atom 
%    and stores the result in HamDG.atomList(:).force
%    
%    See also HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


F = HamDG.fft;
dim = dimDef();

numAtom = length(HamDG.atomList);
numElem = HamDG.numElem;

force = zeros(numAtom, dim);


% ***********************************************************************
% Compute the derivative of the Hartree potential for computing the local
% pseudopotential contribution to Hellmann-Feynman force
% ***********************************************************************
vhartDrv = cell(dim, 1);

totalChargeLocal = cell(HamDG.numElem);
% The contribution of the pseudocCharge is subtracted. So the
% Poisson equation is well defined for neutral system.
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            totalChargeLocal{i, j, k} = ...
                HamDG.density{i, j, k} - HamDG.pseudoCharge{i, j, k};
        end
    end
end

% convert tempVec from 3-dim array into 1-dim vector over global domain
totalCharge = ElemVecToGlobal(totalChargeLocal, HamDG.domain.numGridFine, numElem);

totalChargeFourier = F * totalCharge;

% compute the derivative of the Hartree potential via Fourier transform
gkkFine = F.gkkFine;
idxnz = gkkFine ~= 0;
for d = 1 : dim
    ikFine = F.ikFine{d};

    y = totalChargeFourier;
    y(~idxnz) = 0;
    y(idxnz) = y(idxnz) .* 4*pi ./ gkkFine(idxnz) .* ikFine(idxnz);
    tempVec = real(F' * y);

    % convert tempVec from 1-dim vector into 3-dim array over global domain
    totalChargeLocal = GlobalVecToElem(tempVec, HamDG.domain.numGridFine, numElem);
    vhartDrv{d} = totalChargeLocal;
end


% ***********************************************************************
% Compute the force from local pseudopotential
% ***********************************************************************

% Method 2: Using integration by parts
% This method only uses the value of the local pseudopotential and
% does not use the derivative of the pseudopotential. This is done
% through integration by parts, and the derivative is applied to the
% Coulomb potential evaluated on a uniform grid. 
% 

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            ppMap = HamDG.pseudoListElem{i, j, k};
            for atomIdxcell = keys(ppMap)
                atomIdx = atomIdxcell{1};  % from cell to number
                pseudo = ppMap(atomIdx);
                sp = pseudo.pseudoCharge;
                idx = sp.idx;
                val = sp.val;

                wgt = HamDG.domain.Volume() / HamDG.domain.NumGridTotalFine();
                for d = 1 : dim
                    drv = vhartDrv{d}{i, j, k};
                    res = sum( val(:, 1) .* drv(idx) * wgt );
                    force(atomIdx, d) = force(atomIdx, d) + res;
                end
            end
        end
    end
end


% ***********************************************************************
% compute the force from nonlocal pseudopotential
% ***********************************************************************

% This method only uses the value of the pseudopotential and does not
% use the derivative of the pseudopotential. This is done through
% integration by parts, and the derivative is applied to the basis functions
% evaluated on a LGL grid. See also CalculateDGMatrix.m 

% Loop through the atoms and eigenvecs for the contribution to the force

% loop over atoms and pseudopotentials
numEig = length(HamDG.occupationRate);
for atomIdx = 1 : numAtom
    vnlWeight = HamDG.vnlWeightMap(atomIdx);
    numVnl = length(vnlWeight);
    if numVnl ~= 0    
        resVal = zeros(numEig, numVnl);
        resDrvX = zeros(numEig, numVnl);
        resDrvY = zeros(numEig, numVnl);
        resDrvZ = zeros(numEig, numVnl);
    
        % loop over the elements overlapping with the nonlocal pseudopotenital
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    coefMap = HamDG.vnlCoef{i, j, k};
                    coefDrvXMap = HamDG.vnlDrvCoef{1}{i, j, k};
                    coefDrvYMap = HamDG.vnlDrvCoef{2}{i, j, k};
                    coefDrvZMap = HamDG.vnlDrvCoef{3}{i, j, k};
    
                    localCoef = HamDG.eigvecCoef{i, j, k};
    
                    if isKey(coefMap, atomIdx)
                        coef = coefMap(atomIdx);
                        coefDrvX = coefDrvXMap(atomIdx);
                        coefDrvY = coefDrvYMap(atomIdx);
                        coefDrvZ = coefDrvZMap(atomIdx);
    
                        % skip the calculation if there is no adaptive local
                        % basis function
                        if isempty(coef)
                            continue;
                        end
    
                        % value
                        resVal = resVal + localCoef' * coef;
    
                        % derivative
                        resDrvX = resDrvX + localCoef' * coefDrvX;
                        resDrvY = resDrvY + localCoef' * coefDrvY;
                        resDrvZ = resDrvZ + localCoef' * coefDrvZ;
    
                    end  % found the atom
                end
            end
        end
    
        % Add the contribution to the local force
        % The minus sign comes from integration by parts
        % The 4.0 comes form spin(2.0) and that |l> appears twice (2.0)
        occupationRate = HamDG.occupationRate;  
    
        force(atomIdx, 1) = force(atomIdx, 1) - ...
            4.0 * sum( (occupationRate * vnlWeight) .* resVal .* resDrvX, 'all');
        force(atomIdx, 2) = force(atomIdx, 2) - ...
            4.0 * sum( (occupationRate * vnlWeight) .* resVal .* resDrvY, 'all');
        force(atomIdx, 3) = force(atomIdx, 3) - ...
            4.0 * sum( (occupationRate * vnlWeight) .* resVal .* resDrvZ, 'all');
    end
end
 

% ***********************************************************************
% Give the value to atomList
% ***********************************************************************
for a = 1 : numAtom
    HamDG.atomList(a).force = force(a, :);
end


end