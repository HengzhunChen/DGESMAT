function atomDensity = CalculateAtomDensity(HamDG, ptable)
% HAMILTONIANDG/CALCULATEATOMDENSITY generates atom density.
%
%    atomDensity = CalculateAtomDensity(HamDG, ptable) generates atom
%    density by PeriodTable ptable with respect to domain and atomList of
%    HamDG and stores as a cell structure atomDensity, which is used as 
%    initial guess of density in the SCF iteration.
%
%    See also HamiltonianDG, PeriodTable/CalculateAtomDensity, Atom, 
%    SCFDG/Setup.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


ntotFine = HamDG.domain.NumGridTotalFine();
numAtom = length(HamDG.atomList);
vol = HamDG.domain.Volume();
numElem = HamDG.numElem;
numElemTotal = prod(numElem);

F = HamDG.fft;

% the number of electrons for normalization purpose
nelec = 0;
for i = 1 : numAtom
    atype = HamDG.atomList(i).type;
    if ~ptable.pteMap.isKey(atype)
        error('Cannot find the atom type');
    end
    nelec = nelec + ptable.Zval(atype);
end
if mod(nelec, 2) ~= 0
    error('This is spin-restricted calculation, nelec should be even.');
end

% search for the number of atom types and build a list of atom types
atomTypeSet = [ HamDG.atomList(:).type ];
atomTypeSet = unique(atomTypeSet);

% For each atom type, construct the atomic pseudocharge within the cutoff
% radius starting from the origin in the real space, and construct the
% structure factor. This is done by first generating the density on the
% element-wise grid, convert to the Fourier grid, and then back.

atomDensityR = cell(numElemTotal, 1);

tempAtomDensityG = zeros(ntotFine, 1);

for atype = atomTypeSet
    fakeAtom = Atom();
    fakeAtom.type = atype;
    % NOTE: this should be starting point of the global domain
    fakeAtom.pos = HamDG.domain.posStart;
    
    % compute the atomic density in all elements
    for elemIdx = 1 : numElemTotal
        % HamDG.uniformGridElemFine{elemIdx} starts from the posStart
        % of HamDG.domainElem{elemIdx}
        atomDensityR{elemIdx} = ptable.CalculateAtomDensity(...
                fakeAtom, HamDG.domain, ...
                HamDG.grid.uniformGridElemFine{elemIdx});
    end
    
    % convert atomic density into 1-dim vector
    tempAtomDensityR = ElemVecToGlobal(atomDensityR, HamDG.domain.numGridFine, numElem);
    
    % Multiply with the structure factor on the FFT grid
    % Backward FFT is performed only in the end
    ccvec = zeros(ntotFine, 1);
    ikx = F.ikFine{1};
    iky = F.ikFine{2};
    ikz = F.ikFine{3};

    for i = 1 : numAtom
        if HamDG.atomList(i).type == atype
            xx = HamDG.atomList(i).pos(1);
            yy = HamDG.atomList(i).pos(2);
            zz = HamDG.atomList(i).pos(3);
            
            phase = -( ikx.*xx + iky.*yy + ikz.*zz );
            ccvec = ccvec + exp(phase);
        end
    end
            
    % transfer the atomic charge from real space to Fourier space, and 
    % multiply with the structure factor
    Y = F * tempAtomDensityR;
    
    % make it smoother: AGGREESIVELY truncate components beyond EcutWavefun
    idxnz = (F.gkkFine ./ 2) < HamDG.ecutWavefun;
    tempAtomDensityG = tempAtomDensityG + Y .* ccvec;

    tempAtomDensityG(~idxnz) = 0;
end

X = F' * tempAtomDensityG;
tempAtomDensityR = real(X);

% convert atomic density from 1-dim vector to element-wise
atomDensity = GlobalVecToElem(tempAtomDensityR, HamDG.domain.numGridFine, numElem);

sumrho = sum(tempAtomDensityR);
sumrho = sumrho * vol / ntotFine;

% adjustment should be multiplicative
factor = nelec / sumrho;
for elemIdx = 1 : numElemTotal
    atomDensity{elemIdx} = atomDensity{elemIdx} .* factor;
end

end