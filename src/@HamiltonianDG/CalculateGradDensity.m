function gradDensity = CalculateGradDensity(HamDG)
% HAMILTONIANDG/CALCULATEGRADDENSITY calculates gradient of density.
% 
%    gradDensity = CalculateGradDensity(HamDG) computes gradient of 
%    density which has multi-directions for each dim and is stores as 
%    a cell structure gradDensity. 
% 
%    See also HamiltonianDG, HamiltonianDG/CalculateDensity.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

F = HamDG.fft;

dim = dimDef();
density = HamDG.density;
numGridFine = HamDG.domain.numGridFine;
numElem = HamDG.numElem;
gradDensity = cell(dim, 1);

tempVec = ElemVecToGlobal(density, numGridFine, numElem);

tempVecFFT = F * tempVec;
idxnz = F.gkkFine ~= 0;

for d = 1 : dim
    Y = tempVecFFT;
    ik = F.ikFine{d};
    Y(idxnz) = Y(idxnz) .* ik(idxnz);
    Y(~idxnz) = 0;
    X = real(F' * Y);

    GRho = GlobalVecToElem(X, numGridFine, numElem);
    gradDensity{d} = GRho;
end


end