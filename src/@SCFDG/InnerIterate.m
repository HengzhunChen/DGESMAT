function scfDG = InnerIterate(scfDG, outerIter)
% SCFDG/INNERITERATE Inner SCF iteration of DGDFT.
%
%    scfDG = InnerIterate(scfDG, outerIter) runs the inner iteration of
%    outerIter-th outer SCF iteration and saves the result in SCFDG.
%
%    See also SCFDG, HamiltonianDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = scfDG.numElem;
numElemTotal = prod(numElem);
controlVar = scfDG.controlVar;
DGSolver = scfDG.DGSolver;

isInnerSCFConverged = false;


for innerIter = 1 : controlVar.scfInnerMaxIter
    if isInnerSCFConverged
        break;
    end
    scfDG.scfTotalInnerIter = scfDG.scfTotalInnerIter + 1;
    
    timeIterStart = tic;    
    InfoPrint(0, '\nInner iteration # %d of outer iteration # %d starts. \n\n', ...
        innerIter, outerIter);
    
    % ********************************************************************
    % Update potential and construct/update the DG matrix
    % ********************************************************************
    
    timeStart = tic;
    if innerIter == 1
        % The first inner iteration does not update the potential, and
        % construct the global Hamiltonian matrix from scratch
        
        scfDG.hamDG = CalculateDGMatrix(scfDG.hamDG);
                
    else
        % The consequent inner iterations update the potential in the
        % element, and only update the global Hamiltonian matrix.
        
        % update the potential in the element (and the extended element)
        
        % save the old potential on the LGL grid
        scfDG.vtotLGLSave = scfDG.hamDG.vtotLGL;
        
        % update the local potential on the extended element and on the
        % element
        scfDG = UpdateElemLocalPotential(scfDG);
        
        % save the difference of the potential on the LGL grid into
        % scfDG.vtotLGLSave
        for elemIdx = 1 : numElemTotal
            scfDG.vtotLGLSave{elemIdx} = ...
                scfDG.hamDG.vtotLGL{elemIdx} - scfDG.vtotLGLSave{elemIdx};
        end
                
        % update the DG matrix
        scfDG.hamDG = UpdateDGMatrix(scfDG.hamDG, scfDG.vtotLGLSave);

    end
    timeEnd = toc(timeStart);
    InfoPrint(0, ' Total time for constructing or updating DG matrix is %8f [s] \n\n', timeEnd);
    
    
    % ********************************************************************
    % Evaluate the density matrix
    % 
    % Options: eigs
    % ********************************************************************

    %
    % save the mixing variable first
    %
    
    if scfDG.mix.mixVariable == "density"
        scfDG.mix.mixInnerSave = scfDG.hamDG.density;
    elseif scfDG.mix.mixVariable == "potential"
        scfDG.mix.mixInnerSave = scfDG.hamDG.vtot;
    end
        
    %
    % solve HMat
    % 
    timeSolveHMatStart = tic;
    if DGSolver == "eigs"
        sizeHMat = scfDG.hamDG.sizeHMat;
        numEig = scfDG.hamDG.NumStateTotal();
        if sizeHMat < numEig
            error('number of ALBs is too small and less than number of states.');
        end
        if sizeHMat <= 5000
            % for small matrix, use eig() instead of eigs()
            % convert cells from each element into a matrix
            HMat = scfDG.hamDG.HMat;
            for i = 1 : numElemTotal
                for j = 1 : numElemTotal
                    if isempty(HMat{i, j})
                        numRow = size(HMat{i, i}, 1);
                        numCol = size(HMat{j, j}, 1);
                        HMat{i, j} = zeros(numRow, numCol);
                    end
                end
            end
            HMat = cell2mat(HMat);
            HMat = (HMat + HMat') / 2;                
            [eigVec, eigVal] = eig(HMat);
            eigVal = real(diag(eigVal));
            eigVal = eigVal(1 : numEig);
            [eigVal, id] = sort(eigVal);
            eigVec = eigVec(:, id);
            scfDG.hamDG.eigVal = eigVal;
        else
            % for large matrix, store in sparse mode and use eigs()
            HMat = ElemMatToSparse(scfDG.hamDG);
            HMat = (HMat + HMat') / 2;        
            [eigVec, eigVal] = eigs(HMat, numEig, 'smallestreal');
            scfDG.hamDG.eigVal = diag(eigVal);
        end

        % separate matrix of eigenvectors into cells in each element
        eigVecElem = cell(numElemTotal, 1);
        for elemIdx = 1 : numElemTotal
            basisIdx = scfDG.hamDG.elemBasisIdx{elemIdx};
            eigVecElem{elemIdx} = eigVec(basisIdx, :);
        end
        scfDG.hamDG.eigvecCoef = eigVecElem;    
    else
        error('wrong type of DG_Solver');
    end  
    timeSolveHMatEnd = toc(timeSolveHMatStart);
    InfoPrint(0, ' Total time for solving DG matrix is %8f [s] \n\n', timeSolveHMatEnd);

    
    % ********************************************************************
    % Post processing
    % ********************************************************************
    
    % ----------------------------------------------------------------
    % calculate Harris energy
    % ----------------------------------------------------------------
    
    % compute the occupation rate - specific smearing types dealt
    % with within this function
    [scfDG.hamDG.occupationRate, scfDG.fermi] = ...
        scfDG.CalculateOccupationRate(scfDG.hamDG.eigVal);

    % compute the Harris energy functional
    % NOTE: In computing the Harris energy, the density and the
    % potential must be the INPUT density and potential without ANY
    % update.
    scfDG.EfreeHarris = scfDG.CalculateHarrisEnergy();
    
    % -------------------------------------------------------------
    % calculate the new electron density
    % ------------------------------------------------------------- 
    
    [scfDG.hamDG.density, scfDG.hamDG.densityLGL] = scfDG.hamDG.CalculateDensity();
        
    % -----------------------------------------------------------------
    % Update the output potential, and the KS energy
    % -----------------------------------------------------------------
    
    % Update the Hartree energy and the exchange correlation energy and
    % potential for computing the KS energy.
    % NOTE: vtot should not be updated until finishing the computation of
    % the energies.
    
    if contains(scfDG.XCType, "GGA")
        scfDG.hamDG.gradDensity = CalculateGradDensity(scfDG.hamDG);
    end
    
    [scfDG.Exc, scfDG.hamDG.epsxc, scfDG.hamDG.vxc] = scfDG.hamDG.CalculateXC();
    
    scfDG.hamDG.vhart = scfDG.hamDG.CalculateHartree();    
        
    % Compute the KS energy
    scfDG = scfDG.CalculateKSEnergy();
    
    % Update the total potential AFTER updateing the energy
    
    % No external potential
    
    % Compute the new total potential
    scfDG.hamDG.vtot = scfDG.hamDG.CalculateVtot();    
    
    % --------------------------------------------------------------
    % Compute the a posterior error estimator at every step
    % FIXME This is not used when intra-element parallelization is used
    % --------------------------------------------------------------
    
    if scfDG.userOption.isCalculateAPosterioriEachSCF
        %
        % TODO
        %
        error('Currently option isCalculateAPosterioriEachSCF is not supported !');
    end


    % ********************************************************************
    % Mixing
    % ********************************************************************

    % --------------------------------------------------------------
    % compute the error of the mixing variable
    % --------------------------------------------------------------
    normMixDif = 0;
    normMixOld = 0;
    for elemIdx = 1 : numElemTotal
        if scfDG.mix.mixVariable == "density"
            oldVec = scfDG.mix.mixInnerSave{elemIdx};
            newVec = scfDG.hamDG.density{elemIdx};
            normMixDif = normMixDif + norm(oldVec - newVec)^2;
            normMixOld = normMixOld + norm(oldVec)^2;
        elseif scfDG.mix.mixVariable == "potential"
            oldVec = scfDG.mix.mixInnerSave{elemIdx};
            newVec = scfDG.hamDG.vtot{elemIdx};
            normMixDif = normMixDif + norm(oldVec - newVec)^2;
            normMixOld = normMixOld + norm(oldVec)^2;
        end
    end
    
    normMixDif = sqrt(normMixDif);
    normMixOld = sqrt(normMixOld);
    
    scfDG.scfInnerNorm = normMixDif / normMixOld;
    
    if scfDG.scfInnerNorm < controlVar.scfInnerTolerance
        % converged %
        InfoPrint(0, 'Inner SCF is converged!');
        isInnerSCFConverged = true;
    end
    
    
    % --------------------------------------------------------------
    % Mixing for the inner SCF iteration.
    % --------------------------------------------------------------
        
    if controlVar.scfInnerMaxIter == 1
        % Maximum inner iteration = 1 means there is no distinction of
        % inner/outer SCF. Anderson mixing uses the global history
        numAndersonIter = scfDG.scfTotalInnerIter;
    else
        % If more than one inner iterations is used, then Anderson only
        % uses local history.        
        numAndersonIter = innerIter;
    end
    
    if scfDG.mix.mixVariable == "density"
        if scfDG.mix.mixType == "anderson" || scfDG.mix.mixType == "kerker+anderson"
            [scfDG.hamDG.density, scfDG.mix.dfInnerMat, scfDG.mix.dvInnerMat] = ...
                scfDG.AndersonMix(numAndersonIter, ...
                                  scfDG.mix.mixInnerSave, ...
                                  scfDG.hamDG.density, ...
                                  scfDG.mix.dfInnerMat, ...
                                  scfDG.mix.dvInnerMat);
        else
            error('Invalid mixing type');
        end
    elseif scfDG.mix.mixVariable == "potential"
        if scfDG.mix.mixType == "anderson" || scfDG.mix.mixType == "kerker+anderson"
            [scfDG.hamDG.vtot, scfDG.mix.dfInnerMat, scfDG.mix.dvInnerMat] = ...
                scfDG.AndersonMix(numAndersonIter, ...
                                  scfDG.mix.mixInnerSave, ...
                                  scfDG.hamDG.vtot, ...
                                  scfDG.mix.dfInnerMat, ...
                                  scfDG.mix.dvInnerMat);
        else
            error('Invalid mixing type');
        end
    end
    
    
    %
    % Post processing for the density mixing. Make sure that the density is
    % positive, and compute the potential again. This is only used for
    % density mixing.
    %
    if scfDG.mix.mixVariable == "density"
       sumRho = 0;
       for elemIdx = 1 : numElemTotal
           density = scfDG.hamDG.density{elemIdx};
           
           density = max(density, 0);
           sumRho = sumRho + sum(density, 'all');
       end
       sumRho = sumRho * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
       
       rhofac = scfDG.hamDG.numSpin * scfDG.hamDG.numOccupiedState / sumRho;
              
       % normalize the electron density in the global domain
       for elemIdx = 1 : numElemTotal
           localRho = scfDG.hamDG.density{elemIdx};
           scfDG.hamDG.density{elemIdx} = localRho .* rhofac;
       end
       
       % Update the potential after mixing for the next iteration.
       % This is only used for density mixing
       
       % Compute the exchange-correlation potential and energy from the new
       % density
       hamDG = scfDG.hamDG;
        
       if contains(scfDG.XCType, "GGA")
            hamDG.gradDensity = CalculateGradDensity(hamDG);
       end
       
       [scfDG.Exc, hamDG.epsxc, hamDG.vxc] = CalculateXC(hamDG);
       
       hamDG.vhart = CalculateHartree(hamDG);
       
       % No external potential
       
       % compute the new total  potential
       hamDG.vtot = CalculateVtot(hamDG);
       
       scfDG.hamDG = hamDG;
    end
    
    
    % ********************************************************************
    % Print out the state variables of the current iteration
    % ********************************************************************
    
    scfDG.PrintState();
    
    timeIterEnd = toc(timeIterStart);
    
    InfoPrint(0, ' Time for this inner SCF iteration = %8f [s]\n\n', timeIterEnd);
    
    
end

end