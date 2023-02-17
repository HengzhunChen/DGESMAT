function vhart = CalculateHartree(HamDG)
% HAMILTONIANDG/CALCULATEHARTREE calculates Hartree potential.
%
%    vhart = CalculateHartree(HamDG) computes Hartree potential for each
%    element and stores as a cell vhart.
%
%    See also HamiltonianDG, HamiltonianDG/CalculateVtot.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


F = HamDG.fft;

tempVecLocal = cell(HamDG.numElem);

numElem = HamDG.numElem;
numElemTotal = prod(numElem);

for elemIdx = 1 : numElemTotal
    tempVecLocal{elemIdx} = ...
        HamDG.density{elemIdx} - HamDG.pseudoCharge{elemIdx};
end

% convert tempVec from 3-dim array into 1-dim vector over global domain
tempVec = ElemVecToGlobal(tempVecLocal, HamDG.domain.numGridFine, numElem);

% The contribution of the pseudoCharge is subtracted. So the
% Poisson equation is well defined for neutral system.

y = F * tempVec;
idxnz = F.gkkFine ~= 0;
y(idxnz) = y(idxnz) .* 4 * pi ./ F.gkkFine(idxnz);
y(~idxnz) = 0;
tempVec = real(F' * y);

% convert tempVec from 1-dim vector into 3-dim array over global domain
tempVecLocal = GlobalVecToElem(tempVec, HamDG.domain.numGridFine, numElem);

vhart = tempVecLocal;

end