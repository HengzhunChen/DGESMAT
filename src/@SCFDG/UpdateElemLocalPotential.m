function scfDG = UpdateElemLocalPotential(scfDG)
% SCFDG/UPDATEELEMLOCALPOTENTIAL Update the local potential in the 
%    extended element and the element.
%
%    scfDG = UpdateElemLocalPotential(scfDG) updates
%    (1) scfDG.vecEigSol{:}.hamKS.vtot, i.e., local potential in each
%    extended element;
%    (2) scfDG.hamDG.vtotLGL{:}, i.e., local potential on the LGL grid for
%    each element.
%
%    See also SCFDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = scfDG.numElem;
numElemTotal = prod(numElem);
vtot = scfDG.hamDG.vtot;
vtot = reshape(vtot, numElem);
numGridElem = scfDG.hamDG.grid.numUniformGridElemFine;


% Update of the local potential in each extended element locally.
% The nonlocal potential does not need to be updated.
%
% Also update the local potential on the LGL grid in hamDG
%
% NOTE:
% 1. It is hard coded that the extended element is 1 or (2*bufferSize+1) 
% times the size of the element

bufferSize = scfDG.bufferSize;
bufferIdx = -bufferSize : 1 : bufferSize;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            elemIdx = ElemKeyToIdx(i, j, k, numElem);

            % skip the calculation if there is no adaptive local basis
            % function
            eigSol = scfDG.vecEigSol{elemIdx};
            if eigSol.psi.NumStateTotal() == 0
                continue;
            end
            
            % compute the indices of nearby element
            nearbyIdx = cell(1, 3);
            if numElem(1) >= 3
                tempi = i - 1;
                tempIdx = mod(tempi + bufferIdx, numElem(1));
                nearbyIdx{1} = tempIdx + 1;
            else
                nearbyIdx{1} = i;
            end
            if numElem(2) >= 3
                tempj = j - 1;
                tempIdx = mod(tempj + bufferIdx, numElem(2));
                nearbyIdx{2} = tempIdx + 1;
            else
                nearbyIdx{2} = j;
            end
            if numElem(3) >= 3
                tempk = k - 1;
                tempIdx = mod(tempk + bufferIdx, numElem(3));
                nearbyIdx{3} = tempIdx + 1;
            else
                nearbyIdx{3} = k;
            end
                                                                                                                                                          
            % add vtot to extended element
            vtotExtElemCell = cell(length(nearbyIdx{1}), ...
                length(nearbyIdx{2}), length(nearbyIdx{3}));
            kcell = 0;
            for kk = nearbyIdx{3}
                kcell = kcell + 1;
                jcell = 0;
                for jj = nearbyIdx{2}
                    jcell = jcell + 1;
                    icell = 0;
                    for ii = nearbyIdx{1}
                        icell = icell + 1;
                        vtotTns = reshape(vtot{ii, jj, kk}, numGridElem);
                        vtotExtElemCell{icell, jcell, kcell} = vtotTns;
                    end
                end
            end

            vtotExtElem = cell2mat(vtotExtElemCell);
            vtotExtElem = reshape(vtotExtElem, [], 1);
            
            scfDG.vecEigSol{elemIdx}.hamKS.vtot = vtotExtElem;
        end
    end
end

            
% Modify the potential in the extended element. Current options are 
%
% 1. Add barrier (obsolete)
% 2. Periodize the potential
%
% Numerical results indicate that option 2 seems to be better.

for elemIdx = 1 : numElemTotal
    eigSol = scfDG.vecEigSol{elemIdx};
        
    if scfDG.userOption.isPeriodizePotential
        % get the potential
        vtot = eigSol.hamKS.vtot;
        
        % bring the potential to the vacuum level
        vBubble = scfDG.periodicPotential.vBubble;
        [vBubbleX, vBubbleY, vBubbleZ] = ...
            ndgrid(vBubble{1}, vBubble{2}, vBubble{3});
        vBubbleXYZ = [vBubbleX(:), vBubbleY(:), vBubbleZ(:)];
        eigSol.hamKS.vext = (vtot - 0.0) .* (prod(vBubbleXYZ, 2) - 1.0);
        
        % NOTE:
        % Directly modify the vtot. vext is not used in the
        % matrix-vector multiplication in the eigensolver.
        eigSol.hamKS.vtot = eigSol.hamKS.vtot + eigSol.hamKS.vext;
    end
    
    scfDG.vecEigSol{elemIdx} = eigSol;
end


% Update the potential in element on LGL grid
%
% The local potential on the LGL grid is done by using Fourier
% interpolation from the extended element to the element. Gibbs phenomena
% MAY be there but at least this is better than Lagrange interpolation on a
% uniform grid.
%
% NOTE: The interpolated potential on the LGL grid is taken to be the
% MODIFIED potential with vext on the extended element. Therefore it is
% important that the artificial vext vanishes inside the element. When
% periodization option is used, it can potentially reduce the effect of
% Gibbs phenomena.

for elemIdx = 1 : numElemTotal
    eigSol = scfDG.vecEigSol{elemIdx};            
    vtotExtElem = eigSol.hamKS.vtot;
    numGridExtElemFine = eigSol.fft.domain.numGridFine;
    numLGLGrid = scfDG.hamDG.grid.numLGLGridElem;            

    % interpolate from periodic uniform fine grid to LGL grid
    transferMat = scfDG.PeriodicUniformFineToLGLMat;
    vtotLGLElem = InterpByTransferMat(transferMat, numGridExtElemFine, ...
                                      numLGLGrid, vtotExtElem);
    scfDG.hamDG.vtotLGL{elemIdx} = vtotLGLElem;
end

end