function HamDG = CalculatePseudoPotential(HamDG, ptable)
% HAMILTONIANDG/CALCULATEPSEUDOPOTENTIAL calculates pseudo potential.
%
%    HamDG = CalculatePseudoPotential(HamDG, ptable) computes         
%    (1) number of occupied states and stores as HamDG.numOccupiedState;
%    (2) weight of non-local pseudo potential for each atom in 
%        HamDG.atomList ans stores as HamDG.vnlWeightMap;
%    (3) atomic pseudo potential over each element. For each element, 
%        there is a map containing the list of pseudo potential for those 
%        atoms who overlap with current element. All maps over elements 
%        are stores in cell HamDG.pseudoListElem;
%    (4) local potential generated by atoms for each element, i.e., 
%        computes Gaussian pseudoCharge and save as a cell 
%        HamDG.pseudoCharge, compute short range part for VLocal and 
%        save as HamDG.vLocalSR;
%
%    See also HamiltonianDG, PeriodTable.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numAtom = length(HamDG.atomList);
numElem = HamDG.numElem;
numElemTotal = prod(numElem);

% Very important cleaning, especially when updating the atomic positions
HamDG.pseudoListElem = cell(numElemTotal, 1);


% *********************************************************************
% Atomic information
% *********************************************************************

% Calculate the number of occupied states
nelec = 0;
for i = 1 : numAtom
    atype = HamDG.atomList(i).type;
    if ~ptable.pteMap.isKey(atype)
        error('Cannot find the atom type for atom # %d', i);
    end
    nelec = nelec + ptable.Zval(atype);
end

if mod(nelec, 2) ~= 0
    error('This is a spin-restricted calculation. nelec should be even.');
end

HamDG.numOccupiedState = nelec / HamDG.numSpin;


% *********************************************************************
% Generate the atomic pseudopotentials
% *********************************************************************

% Check that the cutoff radius of the pseudopotential is smaller than the
% length of the element.
minLength = min(HamDG.domainElem{1}.length);
for i = 1 : numAtom
    type = HamDG.atomList(i).type;
    % for the case where there is no nonlocal pseudopotential 
    if ptable.IsNonlocal(type)
        Rzero = ptable.RcutNonlocal(type);
    else
        Rzero = 0;
    end
    
    if Rzero >= minLength
        msg = "In order for the current DG partition to work, " + ...
            "the support of the nonlocal pseudopotential must be smaller than " + ...
            "the length of the element along each dimension. " + ...
            "It is now found for atom " + num2str(i) + ...
            " , which is of type " + num2str(type) + ", Rzero = " + num2str(Rzero) + ...
            " , while the length of element is " + num2str(HamDG.domainElem{1}.length);
        error(msg);
    end
end

% Also prepare the integration weights for constructing the DG matrix
% later.

HamDG.vnlWeightMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
getWeightFlag = zeros(numAtom, 1);

for elemIdx = 1 : numElemTotal
    % for each element, there is a ppMap containing the list of
    % pseudopotential for those atoms who overlap with current element.

    ppMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    gridpos = HamDG.grid.uniformGridElemFine{elemIdx};
    
    for a = 1 : numAtom
        type = HamDG.atomList(a).type;
        Rzero = ptable.RcutVLocalSR(type);
        if ptable.IsNonlocal(type)
            Rzero = max(Rzero, ptable.RcutNonlocal(type));
        end
        
        % compute the minimum distance of this atom to all grid points
        minDist = zeros(dimDef(), 1);
        dmLength = HamDG.domain.length;
        pos = HamDG.atomList(a).pos;
        for d = 1 : dimDef()
            dist = gridpos{d} - pos(d);
            dist = dist - round(dist ./ dmLength(d)) .* dmLength(d);
            minDist(d) = min([abs(dist), Rzero]);
        end

        % If this atom overlaps with this element, compute the
        % pseudopotential
        if norm(minDist) <= Rzero
            pp = struct(...
                'pseudoCharge', [], ...
                'vLocalSR', [], ...
                'vnlList', [] ...
                );
            [pp.vLocalSR, pp.pseudoCharge] = ptable.CalculateVLocal(...
                                HamDG.atomList(a), HamDG.domain, ...
                                HamDG.grid.uniformGridElemFine{elemIdx});
            pp.vnlList = ptable.CalculateNonlocalPP(...
                              HamDG.atomList(a), HamDG.domain, ...
                              HamDG.grid.LGLGridElem{elemIdx});
            ppMap(a) = pp;
            
            if getWeightFlag(a) == 0
                % if atom has not nonlocal pseudopotential part, 
                % set weight = []
                weight = [];
                if ~isempty(pp.vnlList.wgt)
                    weight = pp.vnlList.wgt;
                    getWeightFlag(a) = 1;
                end
                HamDG.vnlWeightMap(a) = weight;
            end
        end
    end  % end for a
    
    HamDG.pseudoListElem{elemIdx} = ppMap;
end


% *********************************************************************
% Local pseudopotential
% *********************************************************************

% compute the Gaussian pseudocharge by summing over contributions from 
% all atoms
HamDG.pseudoCharge = cell(numElemTotal, 1);

sumRho = 0.0;

for elemIdx = 1 : numElemTotal
    ppMap = HamDG.pseudoListElem{elemIdx};
    localVec = zeros(prod(HamDG.grid.numUniformGridElemFine), 1);
    for atomIdxcell = keys(ppMap)
        atomIdx = atomIdxcell{1};
        pp = ppMap(atomIdx);
        sp = pp.pseudoCharge;
        idx = sp.idx;
        val = sp.val;
        localVec(idx) = localVec(idx) + val;
    end
    HamDG.pseudoCharge{elemIdx} = localVec;
    
    sumRho = sumRho + sum(localVec);
end

% compute the sum of the pseudocharge and make adjustment
sumRho = sumRho .* (HamDG.domain.Volume() / HamDG.domain.NumGridTotalFine());

% make adjustments to the pseudocharge
fac = HamDG.numSpin * HamDG.numOccupiedState / sumRho;

for elemIdx = 1 : numElemTotal
    HamDG.pseudoCharge{elemIdx} = HamDG.pseudoCharge{elemIdx} .* fac;
end


% compute the VLocalSR by summing over contributions from all atoms
HamDG.vLocalSR = cell(numElemTotal, 1);
for elemIdx = 1 : numElemTotal
    ppMap = HamDG.pseudoListElem{elemIdx};
    localVec = zeros(prod(HamDG.grid.numUniformGridElemFine), 1);
    for atomIdxCell = keys(ppMap)
        atomIdx = atomIdxCell{1};
        pp = ppMap(atomIdx);
        sp = pp.vLocalSR;
        idx = sp.idx;
        val = sp.val;
        localVec(idx) = localVec(idx) + val;
    end
    HamDG.vLocalSR{elemIdx} = localVec;
end

% ***********************************************************************
% Some other energies and forces
% ***********************************************************************

% [HamDG.Eself, HamDG.EIonSR, HamDG.forceIonSR] = ...
%                         CalculateIonSelfEnergyAndForce(HamDG, ptable);

[HamDG.Eself, HamDG.EIonSR, HamDG.forceIonSR] = ...
    CalculateIonSelfEnergyAndForce(ptable, HamDG.atomList, HamDG.domain);

[HamDG.EVdw, HamDG.forceVdw] = CalculateVdwEnergyAndForce(...
            HamDG.domain, HamDG.atomList, HamDG.VDWType, HamDG.XCType);

HamDG.Eext = 0;

HamDG.forceExt = zeros(length(HamDG.atomList), dimDef());


end