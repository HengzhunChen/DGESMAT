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

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

numElem = scfDG.numElem;
vtot = scfDG.hamDG.vtot;
numGridElem = scfDG.hamDG.grid.numUniformGridElemFine;


% Update of the local potential in each extended element locally.
% The nonlocal potential does not need to be updated.
%
% Also update the local potential on the LGL grid in hamDG
%
% NOTE:
% 1. It is hard coded that the extended element is 1 or 3 times the size of
% the element

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            % skip the calculation if there is no adaptive local basis
            % function
            eigSol = scfDG.vecEigSol{i, j, k};
            if eigSol.psi.NumStateTotal() == 0
                continue;
            end
            
            %  find neighbor of element (i, j, k) to update extended
            %  element (i, j, k)
            vtotExtElemCell = cell(3, 3, 3);
            
            for kk = -1 : 1
                for jj = -1 : 1
                    for ii = -1 : 1
                        % temporarily suppose all 3 dimension have extended element
                        if numElem(1) >= 3
                            idxi = i + ii;
                            if idxi <= 0
                                idxi = idxi + numElem(1);
                            elseif idxi > numElem(1)
                                idxi = mod(idxi, numElem(1));
                            end
                        else
                            idxi = i;
                        end                        
                        if numElem(2) >= 3
                            idxj = j + jj;
                            if idxj <= 0
                                idxj = idxj + numElem(2);
                            elseif idxj > numElem(2)
                                idxj = mod(idxj, numElem(2));
                            end
                        else
                            idxj = j;
                        end
                        if numElem(3) >= 3
                            idxk = k + kk;
                            if idxk <= 0
                                idxk = idxk + numElem(3);
                            elseif idxk > numElem(3)
                                idxk = mod(idxk, numElem(3));
                            end
                        else
                            idxk = k;
                        end
                        
                        vtotTns = reshape(vtot{idxi, idxj, idxk}, numGridElem);
                        vtotExtElemCell{ii+2, jj+2, kk+2} = vtotTns;
                    end
                end
            end
            % remove extra cell
            vtotExtElemTemp = vtotExtElemCell;
            if numElem(1) == 1
                vtotExtElemTemp = vtotExtElemTemp(2, :, :);
            end
            if numElem(2) == 1
                vtotExtElemTemp = vtotExtElemTemp(:, 2, :);
            end
            if numElem(3) == 1
                vtotExtElemTemp = vtotExtElemTemp(:, :, 2);
            end
            
            % update the potential in the extended element
            vtotExtElem = cell2mat(vtotExtElemTemp);
            vtotExtElem = reshape(vtotExtElem, [], 1);
            
            scfDG.vecEigSol{i, j, k}.hamKS.vtot = vtotExtElem;
        end
    end
end

            
% Modify the potential in the extended element. Current options are 
%
% 1. Add barrier
% 2. Periodize the potential
%
% Numerical results indicate that option 2 seems to be better.

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            eigSol = scfDG.vecEigSol{i, j, k};
            
            % Option 1:
            % Add the external barrier potential. CANNOT be used together
            % with periodization option
            if esdfParam.userOption.DG.isPotentialBarrier
                vBarrier = scfDG.potentialBarrier.vBarrier;
                [vBarrierX, vBarrierY, vBarrierZ] = ...
                    ndgrid(vBarrier{1}, vBarrier{2}, vBarrier{3});
                vBarrierXYZ = [vBarrierX(:), vBarrierY(:), vBarrierZ(:)];
                eigSol.hamKS.vext = sum(vBarrierXYZ, 2);
                
                % NOTE:
                % Directly modify the vtot. vext is not used in the
                % matrix-vector multiplication in the eigensolver.
                eigSol.hamKS.vtot = eigSol.hamKS.vtot + eigSol.hamKS.vext;
            end
            
            % Option 2:
            % Periodize the external potential. CANNOT be used together
            % with the barrier potential option
            if esdfParam.userOption.DG.isPeriodizePotential
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
            
            scfDG.vecEigSol{i, j, k} = eigSol;
        end
    end
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

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            eigSol = scfDG.vecEigSol{i, j, k};            
            vtotExtElem = eigSol.hamKS.vtot;
            vtotLGLElem = scfDG.InterpPeriodicUniformFineToLGL(vtotExtElem);
            scfDG.hamDG.vtotLGL{i, j, k} = vtotLGLElem;
        end
    end
end

end