function atomDensity = CalculateAtomDensity(HamKS, ptable)
% HAMILTONIANKS/CALCULATEATOMDENSITY generates atom density.
%
%    atomDensity = CalculateAtomDensity(HamKS, ptable) generates atom
%    density by PeriodTable ptable with respect to domain and atomList of
%    HamKS, atomDensity is used as initial guess of density in the SCF 
%    iteration.
%
%    See also HamiltonianKS, PeriodTable/CalculateAtomDensity, Atom, 
%    SCF/Setup.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


if HamKS.pseudoType == "HGH"
    error('HGH pseudopotential does not yet support the computation of atomic density!');
end

ntotFine = HamKS.domain.NumGridTotalFine();
numAtom = length(HamKS.atomList);
vol = HamKS.domain.Volume();
gridpos = UniformMeshFine(HamKS.domain);

F = HamKS.fft;

% check the number of electrons for normalization purpose
nelec = 0;
for i = 1 : numAtom
    atype = HamKS.atomList(i).type;
    if ~ptable.pteMap.isKey(atype)
        error('Cannot find the atom type');
    end
    nelec = nelec + ptable.Zion(atype);
end
% add extra electron
nelec = nelec + HamKS.numExtraElectron;
if mod(nelec, 2) ~= 0
    error('This is spin-restricted calculation, nelec should be even.');
end

% search for the number of atom types and build a list of atom types 
atomTypeSet = [ HamKS.atomList(:).type ];
atomTypeSet = unique(atomTypeSet);

% for each atom type, construct the atomic pseudocharge within the cutoff
% radius starting from the origin in the real space, and construct the
% structure factor

% Origin-centered atomDensity in the real space and Fourier space
atomDensityG = zeros(ntotFine, 1);

for atype = atomTypeSet
    fakeAtom = Atom();
    fakeAtom.type = atype;
    fakeAtom.pos = HamKS.domain.posStart;
    
    atomDensityR = ptable.CalculateAtomDensity(...
                        fakeAtom, HamKS.domain, gridpos);

    % compute the structure factor
    ccvec = zeros(ntotFine, 1);
    ikx = F.ikFine{1};
    iky = F.ikFine{2};
    ikz = F.ikFine{3};
    
    for i = 1 : numAtom
        if HamKS.atomList(i).type == atype
            xx = HamKS.atomList(i).pos(1);
            yy = HamKS.atomList(i).pos(2);
            zz = HamKS.atomList(i).pos(3);
            
            phase = -( ikx.*xx + iky.*yy + ikz.*zz );
            ccvec = ccvec + exp(phase);
        end
    end
            
    % transfer the atomic charge from real space to Fourier space, and 
    % multiply with the structure factor    
    Y = F * atomDensityR; 
    
    % make it smoother: AGGREESIVELY truncate components beyond EcutWavefunction
    idxnz = (F.gkkFine ./ 2) < HamKS.ecutWavefunction;
    atomDensityG = atomDensityG + Y .* ccvec;
    atomDensityG(~idxnz) = 0;
end
    
% transfer back to the real space and add to HamKS.atomDensity
X = F' * atomDensityG;
atomDensity = real(X);


sumrho = sum(atomDensity);
sumrho = sumrho * vol / ntotFine;
InfoPrint(0, "Sum of atomic density                        = ", sumrho);

% adjustment should be multiplicative
factor = nelec / sumrho;
atomDensity = atomDensity .* factor;

InfoPrint(0, "After adjustment, Sum of atomic density      = ", nelec);

end