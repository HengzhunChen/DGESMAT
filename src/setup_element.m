function VecEigSol = setup_element(esdfParam, ptable)
% SETUP_ELEMENT Setup the element and extended element information
%
%    VecEigSol = setup_element(esdfParam, ptable) returns cells VecEigSol
%    containing EigenSolverKS object over each element with respect to  
%    ESDFInputParam object esdfParam and PeriodTable object ptable.
%    
%    See also EigenSolverKS, HamiltonianKS, Spinor, ESDFInputParam,
%    PeriodTable.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

dm = esdfParam.basic.domain;
numElem = esdfParam.DG.numElem;

VecEigSol = cell(numElem);
VecHamKS  = cell(numElem);
VecPsi    = cell(numElem);


% Element partition
timeStart = tic;
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            
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
                    dmExtElem.length(d)      = dm.length(d) / numElem(d) * 3;
                    dmExtElem.numGrid(d)     = esdfParam.DG.numGridWavefunctionElem(d) * 3;
                    dmExtElem.numGridFine(d) = esdfParam.DG.numGridDensityElem(d) * 3;
                    dmExtElem.posStart(d)    = dm.length(d) / numElem(d) * ( key(d) - 1 );
                else
                    error('numElem(%d) is either 1 or >=3', d);
                end
            end
            
            % Atom
            atomList = esdfParam.basic.atomList;
            numAtom = length(atomList);            
            atomListExtElem = Atom.empty();
            idx = 1;
            for a = 1 : numAtom
                pos = atomList(a).pos;
                if IsInSubdomain(pos, dmExtElem, dm.length)
                    % update the coordinate relative to the extended element
                    pos = pos - floor( (pos - dmExtElem.posStart) ./ dm.length ) ...
                                .* dm.length;
                    atomListExtElem(idx) = Atom(atomList(a).type, ...
                                                pos, ...
                                                atomList(a).vel, ...
                                                atomList(a).force );
                    idx = idx + 1;
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
            numStateTotal = esdfParam.DG.numALBElem(i, j, k) + esdfParam.basic.numUnusedState;
            VecPsi{i, j, k} = Spinor(dmExtElem, numStateTotal, 0.0);
            VecPsi{i, j, k}.wavefun = rand( size(VecPsi{i, j, k}.wavefun) );

            % for test and debug
%             for kk = 1 : numStateTotal
%                 VecPsi{i, j, k}.wavefun(i-1+j-1+k-1+kk, kk) = 1;
%             end
            
            % Hamiltonian KohnSham
            % NOTE: The exchange-correlation type and numExtraState is not used
            % in the extended element calculation
            VecHamKS{i, j, k} = HamiltonianKS(dmExtElem, atomListExtElem, ptable, fftExtElem);
            VecHamKS{i, j, k} = CalculatePseudoPotential(VecHamKS{i, j, k});
            
            % EigenSolverKS
            VecEigSol{i, j, k} = EigenSolverKS(VecHamKS{i, j, k}, VecPsi{i, j, k});
            
        end
    end
end
timeEnd = toc(timeStart);
InfoPrint([0, 1], 'Time for setting up extended element is %8f [s]\n\n', timeEnd);

end