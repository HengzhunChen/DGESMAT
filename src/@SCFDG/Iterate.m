function scfDG = Iterate(scfDG)
% SCFDG/ITERATE Outer SCF iteration of DGDFT.
%
%    See also SCFDG, HamiltonianDG.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

global esdfParam;

numElem = scfDG.numElem;

% **********************************************************************
%                            Preparation                               
% **********************************************************************

dmElem = Domain();
for d = 1 : dimDef()
    dmElem.length(d) = scfDG.domain.length(d) / numElem(d);
    dmElem.numGrid(d) = scfDG.domain.numGrid(d) / numElem(d);
    dmElem.numGridFine(d) = scfDG.domain.numGridFine(d) / numElem(d);
    if numElem(d) > 1
        dmElem.posStart(d) = dmElem.length(d);
    else
        dmElem.posStart(d) = 0;
    end
end

hamDG = scfDG.hamDG;

if contains(scfDG.XCType, "GGA")
    hamDG.gradDensity = CalculateGradDensity(hamDG);
end

[scfDG.Exc, hamDG.epsxc, hamDG.vxc] = CalculateXC(hamDG);

hamDG.vhart = CalculateHartree(hamDG);

% No external potential

hamDG.vtot = CalculateVtot(hamDG);

scfDG.hamDG = hamDG;


% ***********************************************************************
%                            Main Iteration                              
% ***********************************************************************

isSCFConverged = false;
scfDG.scfTotalInnerIter = 0;
controlVar = scfDG.controlVar;

timeTotalStart = tic;

% Total number of SVD basis functions. Determined at the first outer SCF
% and is not changed later. This facilitates the reuse of symbolic
% factorization.
numSVDBasisTotal = 0;

%
% FIXME: Assuming spinor has only one component here
%

for iter = 1 : controlVar.scfOuterMaxIter
    if isSCFConverged && iter >= controlVar.scfOuterMinIter
        break;
    end
    
    % performing each iteration
    PrintBlock(0, ' Outer SCF iteration # ', iter);
    
    timeIterStart = tic;
    
    % *********************************************************************
    % Update the local potential in the extended element and the element.
    %
    % NOTE: The modification of the potential on the extended element
    % to reduce the Gibbs phenomena is now in UpdateElemLocalPotential
    % *********************************************************************
    scfDG = UpdateElemLocalPotential(scfDG);
    
    
    % *********************************************************************
    % Solve the basis functions in the extended element
    % *********************************************************************
    timeBasisStart = tic;
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                eigSol = scfDG.vecEigSol{i, j, k};
                hamDG = scfDG.hamDG;
                
                numLGLGrid = hamDG.grid.numLGLGridElem;
                numGridElemFine = dmElem.numGridFine;                
                
                % skip the interpolation if there is no adaptive local
                % basis function
                if eigSol.psi.NumStateTotal() == 0
                    scfDG.hamDG.basisLGL{i, j, k} = zeros(prod(numLGLGrid), 1);
                    scfDG.hamDG.basisUniformFine{i, j, k} = zeros(prod(numGridElemFine), 1);
                    continue;
                end
                
                % -------------------------------------------------------
                % Solve the basis functions in the extened element
                % -------------------------------------------------------
                
                if controlVar.isPWeigToleranceDynamic
                    % dynamic strategy to control the tolerance
                    if iter == 1
                        % FIXME magic number
                        eigTolNow = 1e-2;
                    else
                        eigTolNow = controlVar.eigTolerance;
                    end
                else
                    % static strategy to control the tolerance
                    eigTolNow = controlVar.eigTolerance;
                end

                eigMaxIter = controlVar.eigMaxIter;
                eigMinTolerance = controlVar.eigMinTolerance;
                
                numEig = eigSol.psi.NumStateTotal() - scfDG.numUnusedState;
                
                % FIXME: multiple choices of solvers for the extended
                % element, should be given in the input file
                if scfDG.PWSolver == "CheFSI"
                    ChebyPW = scfDG.CheFSIPW;
                    % Use CheFSI or LOBPCG on first step
                    if iter <= 1
                        if ChebyPW.firstCycleNum <= 0
                            InfoPrint(0, ' >>>> Calling LOBPCG for ALB generation on extended element ... \n');
                            eigSol = LOBPCGSolve(eigSol, numEig, ...
                                eigMaxIter, eigMinTolerance, eigTolNow);
                        else
                            InfoPrint(0, ' >>>> Calling CheFSI with random guess for ALB generation on extended element ... \n');
                            eigSol = ChebyStepFirst(eigSol, numEig, ...
                                ChebyPW.firstCycleNum, ChebyPW.firstFilterOrder);
                        end
                        InfoPrint(0, '\n');
                    else
                        InfoPrint(0, ' >>>> Calling CheFSI with previous ALBs for generation of new ALBs ... \n');
                        InfoPrint(0, ' >>>> Will carry out %d CheFSI cycles \n', eigMaxIter);
                        for ChebyIter = 1 : eigMaxIter
                            InfoPrint(0, '\n >>>> CheFSI for ALBs : Cycle %d of %d ... \n', ChebyIter, eigMaxIter);
                            eigSol = ChebyStepGeneral(eigSol, numEig,...
                                ChebyPW.generalFilterOrder);
                        end
                        InfoPrint(0, '\n');
                    end
                elseif scfDG.PWSolver == "PPCG"
                    % Use LOBPCG on very first step, i.e., while starting
                    % from random guess
                    if iter <= 1
                        InfoPrint(0, ' >>>> Calling LOBPCG for ALB generation on extended element ... \n');
                        eigSol = LOBPCGSolve(eigSol, numEig, ...
                            eigMaxIter, eigMinTolerance, eigTolNow);
                    else
                        InfoPrint(0, ' >>>> Calling PPCG with previous ALBs for generation of new ALBs ...\n');
                        eigSol = PPCGSolve(eigSol, numEig, eigMaxIter);
                    end
                elseif scfDG.PWSolver == "LOBPCG"
                    % Use LOBPCG to diagonalize
                    eigSol = LOBPCGSolve(eigSol, numEig, ...
                        eigMaxIter, eigMinTolerance, eigTolNow);
                elseif scfDG.PWSolver == "eigs"
                    % Use eigs() function in MATLAB
                    eigSol = eigsSolve(eigSol, numEig, eigMaxIter, eigTolNow);
                else
                    error('Not supported PWSolver type.');
                end
                                
                
                % Assuming that wavefun has only 1 component. This should
                % be changed when spin-polarization is added.
                                
                % compute numBasis in the presence of numUnusedState
                numBasisTotal = eigSol.psi.NumStateTotal() - scfDG.numUnusedState;
                
                elemBasis = zeros(prod(numLGLGrid), numBasisTotal);                
                % NOTE: only one spin component used here
                for ll = 1 : numBasisTotal
                    elemBasis(:, ll) = scfDG.InterpPeriodicUniformToLGL(...
                        eigSol.psi.wavefun(:, ll));
                end

                
                % -----------------------------------------------------------
                % Post processing for the basis functions on the LGL grid.
                % Perform matrix multiplication and threshold the basis 
                % functions for the small matrix.
                %
                % This method might have lower numerical accuracy, but is
                % much more scalable than other known options.
                % ------------------------------------------------------------
                timeStart = tic;
                
                % scale the basis functions by sqrt(weight). This allows
                % the consequent SVD decomposition of the form 
                %
                % X'* W * X
                %
                
                % X(:, i) = sqrt(W) .* X(:, i); 
                sqrtLGLWeight3D = sqrt(scfDG.hamDG.grid.LGLWeight3D);
                sqrtLGLWeight = reshape(sqrtLGLWeight3D, [], 1);
                elemBasis = elemBasis .* repmat(sqrtLGLWeight, 1, numBasisTotal);
                
                MMat = elemBasis' * elemBasis;
                
                [U, S, ~] = svd(MMat, "econ");
                S = diag(S);                
                S = sqrt(S);
                
                % Total number of SVD basis functions. NOTE: Determined at
                % the first outer SCF and is not changed later. This
                % facilitates the reuse of symbolic factorization.
                if iter == 1
                    numSVDBasisTotal = sum( (S ./ S(1)) > controlVar.SVDBasisTolerance);
                else
                    % Reuse the value saved in numSVDBasisTotal
                    InfoPrint(0, 'NOTE: The number of basis functions (after SVD) ');
                    InfoPrint(0, 'is the same as the number in the first SCF iteration.\n');
                end
                
                % Multiply X <- X*U 
                % get the first numSVDBasis which are significant
                U = U(:, 1 : numSVDBasisTotal);
                S = S(1 : numSVDBasisTotal);
                U = U ./ S';
                basis = elemBasis * U;                
                % Unscale the orthogonal basis functions by sqrt of
                % intergration weight
                basis = basis ./ sqrtLGLWeight;
                
                scfDG.hamDG.basisLGL{i, j, k} = basis;
                
                timeEnd = toc(timeStart);
                InfoPrint(0, ' Time for SVD of basis = %8f [s]\n', timeEnd);
                
                scfDG.vecEigSol{i, j, k} = eigSol;                
            end
        end
    end  % end of ALB calculation 
    
    timeBasisEnd = toc(timeBasisStart);
    InfoPrint(0, '\n Time for generating ALB function is %8f [s]\n\n', timeBasisEnd);
    
    
    if scfDG.DGSolver == "CheFSI"
    % -----------------------------------------------------------------
    % Routine for re-orienting eigenvectors based on current basis set
    % -----------------------------------------------------------------
        %
        % TODO 
        %
    end
    
    
    % *********************************************************************
    % Inner SCF iteration 
    %
    % Assemble and diagonalize the DG matrix until convergence is
    % reached for updating the basis functions in the next step.
    % *********************************************************************
    
    % save the mixing variable in the outer SCF iteration
    if scfDG.mix.mixVariable == "density"
        scfDG.mix.mixOuterSave = scfDG.hamDG.density;
    elseif scfDG.mix.mixVariable == "potential"
        scfDG.mix.mixOuterSave = scfDG.hamDG.vtot;
    end
    
    %
    % Main function here
    %
    scfDG = scfDG.InnerIterate(iter);
        
    
    % *********************************************************************
    % Post processing 
    % *********************************************************************
    
    numAtom = length(scfDG.hamDG.atomList);
    scfDG.efreeDifPerAtom = abs(scfDG.Efree - scfDG.EfreeHarris) / numAtom;
    
    % Energy based convergence parameters
    ionDyn = scfDG.ionDyn;
    if iter > 1
        ionDyn.MDscfEbandOld = ionDyn.MDscfEband;
        ionDyn.MDscfEtotOld = ionDyn.MDscfEtot;      
    else
        ionDyn.MDscfEbandOld = 0.0;                
        ionDyn.MDscfEtotOld = 0.0;
    end
    
    ionDyn.MDscfEband = scfDG.Ekin;
    ionDyn.MDscfEbandDiff = abs(ionDyn.MDscfEbandOld - ionDyn.MDscfEband) / numAtom;
    ionDyn.MDscfEtot = scfDG.Etot;
    ionDyn.MDscfEtotDiff = abs(ionDyn.MDscfEtotOld - ionDyn.MDscfEtot) / numAtom;
    scfDG.ionDyn = ionDyn;
    
    % Compute the error of the mixing variable
    normMixDif = 0;
    normMixOld = 0;
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                if scfDG.mix.mixVariable == "density"
                    oldVec = scfDG.mix.mixOuterSave{i, j, k};
                    newVec = scfDG.hamDG.density{i, j, k};
                    normMixDif = normMixDif + norm( oldVec - newVec )^2;
                    normMixOld = normMixOld + norm( oldVec ).^2;
                elseif scfDG.mix.mixVariable == "potential"
                    oldVec = scfDG.mix.mixOuterSave{i, j, k};
                    newVec = scfDG.hamDG.vtot{i, j, k};
                    normMixDif = normMixDif + norm( oldVec - newVec )^2;
                    normMixOld = normMixOld + norm( oldVec )^2;
                end
            end
        end
    end
    
    normMixDif = sqrt(normMixDif);
    normMixOld = sqrt(normMixOld);
    
    scfDG.scfOuterNorm = normMixDif / normMixOld;
    
    InfoPrint([0, 1], 'outer SCF iteration # %d \n', iter);
    InfoPrint([0, 1], "outerSCF: EfreeHarris                 = ", scfDG.EfreeHarris); 
    InfoPrint([0, 1], "outerSCF: Efree                       = ", scfDG.Efree); 
    InfoPrint([0, 1], "outerSCF: norm(out-in)/norm(in) = ", scfDG.scfOuterNorm); 
    InfoPrint([0, 1], "outerSCF: Efree diff per atom   = ", scfDG.efreeDifPerAtom); 

    if scfDG.ionDyn.isUseEnergySCFconvergence
        InfoPrint([0, 1], "outerSCF: MD SCF Etot diff (per atom)           = ", scfDG.ionDyn.MDscfEtotDiff); 
        InfoPrint([0, 1], "outerSCF: MD SCF Eband diff (per atom)          = ", scfDG.ionDyn.MDscfEbandDiff); 
    end

    
    % Check for convergence
    if scfDG.ionDyn.isUseEnergySCFconvergence == 0
        if iter >= 2 ...
                && scfDG.scfOuterNorm < controlVar.scfOuterTolerance ...
                && scfDG.efreeDifPerAtom < controlVar.scfOuterEnergyTolerance
            % converged %
            InfoPrint([0, 1], ' Outer SCF is converged in %d steps ! \n', iter);
            isSCFConverged = true;
        end
    else
        if iter >= 2 ...
                && scfDG.ionDyn.MDscfEtotDiff < scfDG.ionDyn.MDscfEtotDiffTol ...
                && scfDG.ionDyn.MDscfEbandDiff < scfDG.ionDyn.MDscfEbandDiffTol
            % converged via energy criterion
            InfoPrint([0, 1], ' Outer SCF is converged via energy condition in %d step !\n', iter);
            isSCFConverged = true;
        end
    end


    % Potential mixing for the outer SCF iteration. or no mixing at all anymore?
    % It seems that no mixing is the best.
        
    timeIterEnd = toc(timeIterStart);
    InfoPrint(0, ' Time for this SCF iteration = %8f [s]\n', timeIterEnd);
    
end  % -----------------  end of for iter -------------------------------


timeTotalEnd = toc(timeTotalStart);

InfoPrint([0, 1], '\nTotal time for all SCF iterations = %8f [s]\n', timeTotalEnd);
if scfDG.ionDyn.ionDynIter >= 1
    InfoPrint(0, ' Ion dynamics iteration %d : ', scfDG.ionDyn.scfdg_ion_dyn_iter);
end

if isSCFConverged
    InfoPrint(0, ' Total number of outer SCF steps for SCF convergence = %d \n', iter-1);
else
    InfoPrint(0, ' Total number of outer SCF steps (SCF not converged) = %d \n', controlVar.scfOuterMaxIter);
end


% *********************************************************************
% Calculate the VDW contribution and the force
% *********************************************************************

timeForceStart = tic;
if scfDG.DGSolver == "eigs"
    InfoPrint(0, '\n Computing forces using eigenvectors ...\n');
    scfDG.hamDG = scfDG.hamDG.CalculateForce();
elseif scfDG.DGSolver == "CheFSI"
    if ~scfDG.CheFSIDG.isUseCompSubspace
        InfoPrint(0, '\n Computing forces using eigenvectors ...\n');
        scfDG.hamDG = scfDG.hamDG.CalculateForce();
    else
        %
        % TODO
        %
    end
elseif scfDG.DGSolver == "PEXSI"
    % 
    % TODO
    %
end
timeForceEnd = toc(timeForceStart);
InfoPrint(0, ' Time for computing the force is %8f [s] \n\n', timeForceEnd);


% Calculate the VDW energy
if scfDG.VDWType == "DFT-D2"
    [scfDG.Evdw, scfDG.forceVdw] = CalculateVDW(scfDG.hamDG, scfDG.VDWType);
    % update energy
    scfDG.Etot = scfDG.Etot + scfDG.Evdw;
    scfDG.Efree = scfDG.Efree + scfDG.Evdw;
    scfDG.EfreeHarris = scfDG.EfreeHarris + scfDG.Evdw;
    scfDG.Ecor = scfDG.Ecor + scfDG.Evdw;
    
    % update force
    for i = 1 : length(scfDG.hamDG.atomList)
        scfDG.hamDG.atomList(i).force = scfDG.hamDG.atomList(i).force + scfDG.forceVdw(i, :);
    end
end


% *********************************************************************
% Output information after SCF
% *********************************************************************

%
% Print out the Energy
%
PrintBlock([0, 1], 'Energy');

InfoPrint([0, 1], 'NOTE:  Ecor  = Exc - EVxc - Ehart -Eself + EVdw \n');
InfoPrint([0, 1], '       Etot  = Ekin + Ecor \n');
InfoPrint([0, 1], '       Efree = Etot + Entropy \n \n');

InfoPrint([0, 1], "! EfreeHarris       = ",  scfDG.EfreeHarris, "[au]");
InfoPrint([0, 1], "! Etot              = ",  scfDG.Etot, "[au]");
InfoPrint([0, 1], "! Efree             = ",  scfDG.Efree, "[au]");
InfoPrint([0, 1], "! EVdw              = ",  scfDG.Evdw, "[au]"); 
InfoPrint([0, 1], "! Fermi             = ",  scfDG.fermi, "[au]");

InfoPrint([0, 1], '\n  Convergence information : \n');
InfoPrint([0, 1], "! norm(out-in)/norm(in) = ",  scfDG.scfOuterNorm); 
InfoPrint([0, 1], "! Efree diff per atom   = ",  scfDG.efreeDifPerAtom, "[au]"); 

if scfDG.ionDyn.isUseEnergySCFconvergence
    InfoPrint([0, 1], "! MD SCF Etot diff (per atom)  = ",  scfDG.ionDyn.MDscfEtotDiff, "[au]"); 
    InfoPrint([0, 1], "! MD SCF Eband diff (per atom) = ",  scfDG.ionDyn.MDscfEbandDiff, "[au]"); 
end    
   
%
% Print out the force
%
PrintBlock([0, 1], 'Atomic Force');
numAtom = length(scfDG.hamDG.atomList);
forceCM = zeros(1, 3);
for i = 1 : numAtom
    InfoPrint([0, 1], 'atom', i, 'force', scfDG.hamDG.atomList(i).force);
    forceCM = forceCM + scfDG.hamDG.atomList(i).force;
end

InfoPrint([0, 1], '\n');
InfoPrint([0, 1], 'force for centroid  : ', forceCM);
InfoPrint([0, 1], 'Max force magnitude : ', MaxForce(scfDG.hamDG.atomList));
InfoPrint([0, 1], '\n');


%
% Output the atomic structure and other information for describing density, 
% basis function, etc as output file.
%

% output structure information as file
if esdfParam.userOption.general.isOutputStructInfo
    domain = scfDG.domain;
    atomList = scfDG.hamDG.atomList;
    save('STRUCTURE.mat', 'domain', 'atomList');
end

% output restarting information as file

% output the density on the uniform fine grid
if esdfParam.userOption.general.isOutputDensity
    uniformGridElemFine = scfDG.hamDG.grid.uniformGridElemFine;  % cell
    density = scfDG.hamDG.density;  % cell
    save(scfDG.restart.DensityFileName, 'uniformGridElemFine', 'density');
end

% output wavefunction in the extended element
if esdfParam.userOption.DG.isOutputWfnExtElem
    gridPosExtElem = cell(numElem);
    wfnExtElem = cell(numElem);
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                dm = scfDG.vecEigSol{i, j, k}.hamKS.domain;
                gridPosExtElem{i, j, k} = UniformMesh(dm);
                wfnExtElem{i, j, k} = scfDG.vecEigSol{i, j, k}.psi.wavefun;
            end
        end
    end
    save(scfDG.restart.WfnExtElemFileName, 'gridPosExtElem', 'wfnExtElem');
end


end