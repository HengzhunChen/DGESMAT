function scfDG = InnerIterate(scfDG, outerIter)
% SCFDG/INNERITERATE Inner SCF iteration of DGDFT.
%
%    scfDG = InnerIterate(scfDG, outerIter) runs the inner iteration of
%    outerIter-th outer SCF iteration and saves the result in SCFDG.
%
%    See also SCFDG, HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

numElem = scfDG.numElem;
controlVar = scfDG.controlVar;
DGSolver = scfDG.DGSolver;
CheFSIDG = scfDG.CheFSIDG;

isInnerSCFConverged = false;


for innerIter = 1 : controlVar.scfInnerMaxIter
    if isInnerSCFConverged
        break;
    end
    scfDG.scfTotalInnerIter = scfDG.scfTotalInnerIter + 1;
    
    timeIterStart = tic;    
    InfoPrint(0, '\n Inner iteration # %d of outer iteration # %d starts. \n\n', innerIter, outerIter);
    
    % ********************************************************************
    % Update potential and construct/update the DG matrix
    % ********************************************************************
    
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
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    scfDG.vtotLGLSave{i, j, k} = ...
                        scfDG.hamDG.vtotLGL{i, j, k} - scfDG.vtotLGLSave{i, j, k};
                end
            end
        end
                
        % update the DG matrix
        scfDG.hamDG = UpdateDGMatrix(scfDG.hamDG, scfDG.vtotLGLSave);

    end
    
    
    % ********************************************************************
    % Evaluate the density matrix
    % 
    % Options: eigs, PEXSI(TODO), CheFSI(TODO)
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
    if DGSolver == "CheFSI"
        timeStart = tic;

        if scfDG.ionDyn.ionDynIter ~= 0
            if CheFSIDG.isUseCompSubspace
                %
                % TODO
                %
            else
                %
                % TODO
                %
            end
        else
            % 0th MD / Geometry Optimization Step (or static calculation)
            if outerIter == 1
                %
                % TODO
                %
            elseif (outerIter > 1 && outerIter <= CheFSIDG.secondOuterIter)
                %
                % TODO
                %
            else
                if CheFSIDG.isUseCompSubspace
                    %
                    % TODO
                    %
                else
                    %
                    % TODO
                    %
                end
            end
        end

        timeEnd = toc(timeStart);
        if CheFSI.isUseCompSubspace
            InfoPrint(0, '\n Total time for Complementary Subspace method is %8f [s] \n\n', timeEnd);
        else
            InfoPrint(0, '\n Total time for diag DG matrix via Chebyshev filtering is %8f [s] \n\n', timeEnd);
        end

        % save data back to scfDG
        %
        % TODO
        %

    elseif DGSolver == "eigs"                 
        % convert cells from each element into a matrix
        HMat = scfDG.hamDG.ElemMatToSparse();
        
        % TODO: used last result as initial value ?
        numEig = scfDG.hamDG.NumStateTotal();
        HMat = (HMat + HMat') / 2;        
        [eigVec, eigVal] = eigs(HMat, numEig, 'smallestreal');

        scfDG.hamDG.eigVal = diag(eigVal);

        % separate matrix of eigenvectors into cells in each element
        eigVecElem = cell(numElem);
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    basisIdx = scfDG.hamDG.elemBasisIdx{i, j, k};
                    eigVecElem{i, j, k} = eigVec(basisIdx, :);
                end
            end
        end
        scfDG.hamDG.eigvecCoef = eigVecElem;
    
    elseif DGSolver == "PEXSI"
        %
        % TODO
        %
    else
        error('wrong type of DG_Solver')
    end  
        
    
    % ********************************************************************
    % Post processing
    % ********************************************************************

    if DGSolver ~= "PEXSI"

    scfDG.Evdw = 0;
    
    % ----------------------------------------------------------------
    % calculate Harris energy
    % ----------------------------------------------------------------
    
    if DGSolver == "eigs"
        % compute the occupation rate - specific smearing types dealt
        % with within this function
        [scfDG.hamDG.occupationRate, scfDG.fermi] = ...
            scfDG.CalculateOccupationRate(scfDG.hamDG.eigVal);
    
        % compute the Harris energy functional
        % NOTE: In computing the Harris energy, the density and the
        % potential must be the INPUT density and potential without ANY
        % update.
        scfDG.EfreeHarris = scfDG.CalculateHarrisEnergy();
    elseif DGSolver == "CheFSI"
        if CheFSIDG.isUseCompSubspace
            % calculate Harris energy without computing the occupations
            scfDG.EfreeHarris = scfDG.CalculateHarrisEnergy();
        else
            % compute the occupation rate - specific smearing types dealt
            % with within this function
            [scfDG.hamDG.occupationRate, scfDG.fermi] = ...
                scfDG.CalculateOccupationRate(scfDG.hamDG.eigVal);
        
            % compute the Harris energy functional
            % NOTE: In computing the Harris energy, the density and the
            % potential must be the INPUT density and potential without ANY
            % update.
            scfDG.EfreeHarris = scfDG.CalculateHarrisEnergy();
        end
    end
    
    % -------------------------------------------------------------
    % calculate the new electron density
    % ------------------------------------------------------------- 
    
    if DGSolver == "eigs"
        [scfDG.hamDG.density, scfDG.hamDG.densityLGL] = scfDG.hamDG.CalculateDensity();
    elseif DGSolver == "CheFSI"
        if CheFSIDG.isUseCompSubspace
            % density calculation for complementary subspace method
            
            %
            % TODO
            %       
        else
            avgNumALB = scfDG.hamDG.NumBasisTotal() / prod(numElem);  % average no. of ALBs per element
            numStateTotal = scfDG.hamDG.NumStateTotal();
            if avgNumALB < numStateTotal
                %
                % TODO
                %
            else
                [scfDG.hamDG.density, scfDG.hamDG.densityLGL] = scfDG.hamDG.CalculateDensity();
            end
        end
    end
        
    % -----------------------------------------------------------------
    % Update the output potential, and the KS and second order accurate
    % energy
    % -----------------------------------------------------------------
    
    % Update the Hartree energy and the exchange correlation energy and
    % potential for computing the KS energy and the second order energy.
    % NOTE: vtot should not be updated until finishing the computation of
    % the energies.
    
    if contains(scfDG.XCType, "GGA")
        scfDG.hamDG.gradDensity = CalculateGradDensity(scfDG.hamDG);
    end
    
    [scfDG.Exc, scfDG.hamDG.epsxc, scfDG.hamDG.vxc] = scfDG.hamDG.CalculateXC();
    
    scfDG.hamDG.vhart = scfDG.hamDG.CalculateHartree();    
        
    % Compute the second order accurate energy functional.
    % NOTE: In computing the second order energy, the density and the
    % potential must be the OUTPUT density and potential without ANY
    % MIXING.
    scfDG.EfreeSecondOrder = scfDG.CalculateSecondOrderEnergy();
    
    % Compute the KS energy
    scfDG = scfDG.CalculateKSEnergy();
    
    % Update the total potential AFTER updateing the energy
    
    % No external potential
    
    % Compute the new total potential
    scfDG.hamDG.vtot = scfDG.hamDG.CalculateVtot();
    
    
    % -----------------------------------------------------------
    % compute the force at every step
    % -----------------------------------------------------------
    if esdfParam.userOption.DG.isCalculateForceEachSCF
        %
        % TODO
        %
    end
    
    
    % --------------------------------------------------------------
    % Compute the a posterior error estimator at every step
    % FIXME This is not used when intra-element parallelization is used
    % --------------------------------------------------------------
    
    if esdfParam.userOption.DG.isCalculateAPosterioriEachSCF
        %
        % TODO
        %
    end

    end 
    % end of post processing for diagonalization method (i.e. DGSolver ~= PEXSI)
    

    % ********************************************************************
    % Mixing
    % ********************************************************************

    % --------------------------------------------------------------
    % compute the error of the mixing variable
    % --------------------------------------------------------------
    normMixDif = 0;
    normMixOld = 0;
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                if scfDG.mix.mixVariable == "density"
                    oldVec = scfDG.mix.mixInnerSave{i, j, k};
                    newVec = scfDG.hamDG.density{i, j, k};
                    
                    normMixDif = normMixDif + norm(oldVec - newVec)^2;
                    normMixOld = normMixOld + norm(oldVec)^2;
                elseif scfDG.mix.mixVariable == "potential"
                    oldVec = scfDG.mix.mixInnerSave{i, j, k};
                    newVec = scfDG.hamDG.vtot{i, j, k};
                    
                    normMixDif = normMixDif + norm(oldVec - newVec)^2;
                    normMixOld = normMixOld + norm(oldVec)^2;
                end
            end
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
       for k = 1 : numElem(3)
           for j = 1 : numElem(2)
               for i = 1 : numElem(1)
                   density = scfDG.hamDG.density{i, j, k};
                   
                   density = max(density, 0);
                   sumRho = sumRho + sum(density, 'all');
               end
           end
       end
       sumRho = sumRho * scfDG.domain.Volume() / scfDG.domain.NumGridTotalFine();
       
       rhofac = scfDG.hamDG.numSpin * scfDG.hamDG.numOccupiedState / sumRho;
              
       % normalize the electron density in the global domain
       for k = 1 : numElem(3)
           for j = 1 : numElem(2)
               for i = 1 : numElem(1)
                   localRho = scfDG.hamDG.density{i, j, k};
                   scfDG.hamDG.density{i, j, k} = localRho .* rhofac;
               end
           end
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