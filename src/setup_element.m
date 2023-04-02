function VecEigSol = setup_element(esdfParam, ptable)
% SETUP_ELEMENT Setup the element and extended element information
%
%    VecEigSol = setup_element(esdfParam, ptable) returns a cell VecEigSol
%    containing EigenSolverKS object over each element with respect to  
%    ESDFInputParam object esdfParam and PeriodTable object ptable.
%    
%    See also EigenSolverKS, HamiltonianKS, Spinor, ESDFInputParam,
%    PeriodTable.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = esdfParam.DG.numElem;
numElemTotal = prod(numElem);

VecEigSol = cell(numElemTotal, 1);

% Element partition
timeStart = tic;
bufferSize = esdfParam.DG.bufferSize;
extendSize = 2 * bufferSize + 1;

% NOTE: parfor is used here
parfor elemIdx = 1 : numElemTotal
    [i, j, k] = ElemIdxToKey(elemIdx, numElem);

    dm = esdfParam.basic.domain;
            
    % Setup domain in the extended element
    dmExtElem = Domain();
    key = [i-1, j-1, k-1];
    for d = 1 : dimDef()
        % assume the global domain starts from 0.0
        if numElem(d) == 1 
            dmExtElem.length(d)      = dm.length(d);
            dmExtElem.numGrid(d)     = esdfParam.DG.numGridWavefunctionElem(d);
            dmExtElem.numGridFine(d) = esdfParam.DG.numGridDensityElem(d);
            dmExtElem.posStart(d)    = 0.0;
        elseif numElem(d) >= 3
            dmExtElem.length(d)      = dm.length(d) / numElem(d) * extendSize;
            dmExtElem.numGrid(d)     = esdfParam.DG.numGridWavefunctionElem(d) * extendSize;
            dmExtElem.numGridFine(d) = esdfParam.DG.numGridDensityElem(d) * extendSize;
            dmExtElem.posStart(d)    = dm.length(d) / numElem(d) * ( key(d) - bufferSize );
        else
            error('numElem(%d) is either 1 or >=3', d);
        end
    end
    
    % Atom
    atomList = esdfParam.basic.atomList;
    numAtom = length(atomList);            
    atomListExtElem = Atom.empty();
    idx = 0;
    posList = nan(numAtom, 3);
    for a = 1 : numAtom
        pos = atomList(a).pos;
        if IsInSubdomain(pos, dmExtElem, dm.length)
            % update the coordinate relative to the extended element
            pos = pos - floor( (pos - dmExtElem.posStart) ./ dm.length ) ...
                        .* dm.length;

            % check for multi-defined atoms under periodic condition of 
            % the extended element
            % NOTE: this MAY be removed when other boundary condition is used
            repeatFlag = checkRepeatAtom(pos, posList(1:idx, :), dmExtElem, numElem);
            if repeatFlag
                continue;
            end
            % if atom not locates on the boundary, add into atomList
            idx = idx + 1;
            posList(idx, :) = pos;    
            atomListExtElem(idx) = Atom(atomList(a).type, ...
                                        pos, ...
                                        atomList(a).vel, ...
                                        atomList(a).force );
        end  % atom add in the extended element
    end
    
    % Fourier
    % All extended elements share the same Fourier structure
    % NOTE: Fourier fftExtElem can only initialized once and each
    % element use the same structure, here we initilize over each 
    % element for simplicity and for the fact that Domain object 
    % in each fftExtElem will have different start positions 

    fftExtElem = Fourier();
    fftExtElem = Initialize(fftExtElem, dmExtElem);
    fftExtElem = InitializeFine(fftExtElem, dmExtElem);
    
    % Wavefunction
    % ONLY restricted case is considered here
    numStateTotal = esdfParam.DG.numALBElem(elemIdx) + esdfParam.basic.numUnusedState;
    psi = Spinor(dmExtElem, numStateTotal, 0.0);
    psi.wavefun = rand( size(psi.wavefun) );

    % for test and debug
%     for kk = 1 : numStateTotal
%         psi.wavefun(i-1+j-1+k-1+kk, kk) = 1;
%     end
    
    % Hamiltonian KohnSham
    % NOTE: The exchange-correlation type and numExtraState is not used
    % in the extended element calculation
    hamKS = HamiltonianKS(esdfParam, dmExtElem, atomListExtElem, ptable, fftExtElem);
    hamKS = CalculatePseudoPotential(hamKS);
    
    % EigenSolverKS
    VecEigSol{elemIdx} = EigenSolverKS(hamKS, psi);
end

% Check whether number of ALBs in each element is enough
% NOTE: this is important for SCF convergence and need further
% consideration with numALB and numState each element for 
% large scale systems in future. 
for elemIdx = 1 : numElemTotal
    numALB = esdfParam.DG.numALBElem(elemIdx);
    numOccupiedState = VecEigSol{elemIdx}.hamKS.numOccupiedState;

    if numALB < numOccupiedState
        [i, j, k] = ElemIdxToKey(elemIdx, numElem);
        msg = "number of ALBs=" + num2str(numALB) + " in element " + ...
            "(" + num2str(i) + ", " + num2str(j) + ", " + num2str(k) + ")" + ...
            " is less than number of OccupiedStates=" + num2str(numOccupiedState) + ...
            " over corresponding extended element." + ...
            " This MAY have significant impact to SCF convergence.";
        warning(msg);
    end
end

timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for setting up extended element is %8f [s]\n\n', timeEnd);

end


function repeatFlag = checkRepeatAtom(pos, posList, dmExtElem, numElem)
% CHECKREPEATATOM subroutine used to check whether an atom at position pos
%    or its periodic images have been defined in the atom list posList 
%    under periodic condition of extended element dmExtElem.

    % check whether the atom on the boundary
    posImage = cell(3, 1);
    for d = 1 : 3
        if numElem(d) >= 3
            if abs(pos(d) - dmExtElem.posStart(d)) < 1e-14
                posImage{d} = [pos(d), pos(d) + dmExtElem.length(d)];
            elseif abs(pos(d) - dmExtElem.posStart(d) - dmExtElem.length(d)) < 1e-14
                posImage{d} = [dmExtElem.posStart(d), pos(d)];
            else
                posImage{d} = pos(d);
            end
        else
            posImage{d} = pos(d);
        end
    end

    % generate position of periodic image 
    [posImageX, posImageY, posImageZ] = ndgrid(posImage{1}, posImage{2}, posImage{3});
    posImages = [posImageX(:), posImageY(:), posImageZ(:)];
    
    % if atom locates at the boudary, check whether its periodic images 
    % already in the atomList
    repeatFlag = false;
    if size(posImages, 1) ~= 1
        for ii = 1 : size(posImages, 1)
            tempPos = posImages(ii, :);
            dist = sum(abs(tempPos - posList), 2);
            if find(dist < 1e-14)
                repeatFlag = true;
            end
        end
    end

end
