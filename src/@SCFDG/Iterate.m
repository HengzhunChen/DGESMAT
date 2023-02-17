function scfDG = Iterate(scfDG)
% SCFDG/ITERATE Outer SCF iteration of DGDFT.
%
%    See also SCFDG, HamiltonianDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


numElem = scfDG.numElem;
numElemTotal = prod(numElem);

% **********************************************************************
%                            Preparation                               
% **********************************************************************

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
numSVDBasisTotal = zeros(numElemTotal, 1);

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
    
    hamDG = scfDG.hamDG;
    basisLGLNew = cell(numElemTotal, 1);

    vecEigSol = scfDG.vecEigSol;
    vecEigSolNew = cell(numElemTotal, 1);

    numLGLGrid = hamDG.grid.numLGLGridElem;
    numUnusedState = scfDG.numUnusedState;

    PWSolver = scfDG.PWSolver;
    ChebyPW = scfDG.CheFSIPW;
    PPCGsbSize = scfDG.PPCGsbSize;

    PeriodicUniformToLGLMat = scfDG.PeriodicUniformToLGLMat;
    numGridExtElem = vecEigSol{1}.fft.domain.numGrid;

    sqrtLGLWeight3D = sqrt(scfDG.hamDG.grid.LGLWeight3D);
    sqrtLGLWeight = reshape(sqrtLGLWeight3D, [], 1);

    % NOTE: parfor is used here
    parfor elemIdx = 1 : numElemTotal
        eigSol = vecEigSol{elemIdx};
                
        % skip the interpolation if there is no adaptive local
        % basis function
        if eigSol.psi.NumStateTotal() == 0
            basisLGLNew{elemIdx} = zeros(prod(numLGLGrid), 1);
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
        
        numEig = eigSol.psi.NumStateTotal() - numUnusedState;
        
        % FIXME: multiple choices of solvers for the extended element, 
        % should be given in the input file
        if PWSolver == "CheFSI"
            % Use CheFSI or LOBPCG on first step
            if iter <= 1
                if ChebyPW.firstCycleNum <= 0
                    InfoPrint(0, ' >>>> Calling LOBPCG with random guess for ALB generation \n');
                    eigSol = LOBPCGSolve(eigSol, numEig, ...
                        eigMaxIter, eigMinTolerance, eigTolNow);
                else
                    InfoPrint(0, ' >>>> Calling CheFSI with random guess for ALB generation \n');
                    eigSol = ChebyStepFirst(eigSol, numEig, ...
                        ChebyPW.firstCycleNum, ChebyPW.firstFilterOrder);
                end
            else
                InfoPrint(0, ' >>>> Calling CheFSI with previous ALBs for generation of new ALBs \n');
                InfoPrint(0, ' >>>> Will carry out %d CheFSI cycles \n', eigMaxIter);
                for ChebyIter = 1 : eigMaxIter
                    InfoPrint(0, '\n CheFSI for ALBs : Cycle %d of %d ... \n', ChebyIter, eigMaxIter);
                    eigSol = ChebyStepGeneral(eigSol, numEig, ...
                        ChebyPW.generalFilterOrder);
                end
            end
        elseif PWSolver == "PPCG"
            % Use LOBPCG on very first step, i.e., while starting from 
            % random guess
            if iter <= 1
                % Call LOBPCG for ALB generation on extended element
                eigSol = LOBPCGSolve(eigSol, numEig, eigMaxIter, ...
                    eigMinTolerance, eigTolNow);
            else
                % Call PPCG with previous ALBs for generation of new ALBs
                eigSol = PPCGSolve(eigSol, numEig, eigMaxIter, PPCGsbSize);
            end
        elseif PWSolver == "LOBPCG"
            % Use LOBPCG to diagonalize
            eigSol = LOBPCGSolve(eigSol, numEig, eigMaxIter, ...
                eigMinTolerance, eigTolNow);
        elseif PWSolver == "eigs"
            % Use eigs() function in MATLAB
            eigSol = eigsSolve(eigSol, numEig, eigMaxIter, eigTolNow);
        else
            error('Not supported PWSolver type.');
        end
                        

        % Assuming that wavefun has only 1 component. This should
        % be changed when spin-polarization is added.
                        
        % compute numBasis in the presence of numUnusedState
        numBasisTotal = eigSol.psi.NumStateTotal() - numUnusedState;
        
        % interpolate from periodic uniform grid to LGL grid
        elemBasis = zeros(prod(numLGLGrid), numBasisTotal);                
        % NOTE: only one spin component used here
        for ll = 1 : numBasisTotal
            elemBasis(:, ll) = InterpByTransferMat(PeriodicUniformToLGLMat, ...
                numGridExtElem, numLGLGrid, eigSol.psi.wavefun(:, ll));
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
        elemBasis = elemBasis .* repmat(sqrtLGLWeight, 1, numBasisTotal);
        
        MMat = elemBasis' * elemBasis;
        [U, S, ~] = svd(MMat, "econ");
        S = diag(S);                
        S = sqrt(S);
        
        % Total number of SVD basis functions. NOTE: Determined at the 
        % first outer SCF and is not changed later. This facilitates the 
        % reuse of symbolic factorization.
        if iter == 1
            numSVDBasis = sum( (S ./ S(1)) > controlVar.SVDBasisTolerance);
            numSVDBasisTotal(elemIdx) = numSVDBasis;            
        else
            % Reuse the value saved in numSVDBasisTotal
            numSVDBasis = numSVDBasisTotal(elemIdx);
            InfoPrint(0, 'NOTE: The number of basis functions (after SVD) ');
            InfoPrint(0, 'is the same as the number in the first SCF iteration.\n');
        end

        % Multiply X <- X*U 
        % get the first numSVDBasis which are significant
        U = U(:, 1 : numSVDBasis);
        S = S(1 : numSVDBasis);
        U = U ./ S';
        basis = elemBasis * U;                
        % Unscale the orthogonal basis functions by sqrt of integration weight
        basis = basis ./ sqrtLGLWeight;
        
        basisLGLNew{elemIdx} = basis;
        vecEigSolNew{elemIdx} = eigSol;                
        
        timeEnd = toc(timeStart);
        InfoPrint(0, ' Time for SVD of basis = %8f [s]\n', timeEnd);

    end  % end of ALB calculation

    scfDG.hamDG.basisLGL = basisLGLNew;
    scfDG.vecEigSol = vecEigSolNew;
    
    timeBasisEnd = toc(timeBasisStart);
    InfoPrint(0, '\n Time for generating ALB function is %8f [s]\n\n', timeBasisEnd);
        
    
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
        
    % Compute the error of the mixing variable
    normMixDif = 0;
    normMixOld = 0;
    for elemIdx = 1 : numElemTotal
        if scfDG.mix.mixVariable == "density"
            oldVec = scfDG.mix.mixOuterSave{elemIdx};
            newVec = scfDG.hamDG.density{elemIdx};
            normMixDif = normMixDif + norm( oldVec - newVec )^2;
            normMixOld = normMixOld + norm( oldVec ).^2;
        elseif scfDG.mix.mixVariable == "potential"
            oldVec = scfDG.mix.mixOuterSave{elemIdx};
            newVec = scfDG.hamDG.vtot{elemIdx};
            normMixDif = normMixDif + norm( oldVec - newVec )^2;
            normMixOld = normMixOld + norm( oldVec )^2;
        end
    end
    
    normMixDif = sqrt(normMixDif);
    normMixOld = sqrt(normMixOld);
    
    scfDG.scfOuterNorm = normMixDif / normMixOld;
    
    InfoPrint([0, 1], 'outer SCF iteration # %d \n', iter);
    InfoPrint([0, 1], "outerSCF: EfreeHarris                 = ", scfDG.EfreeHarris); 
    InfoPrint([0, 1], "outerSCF: Efree                       = ", scfDG.Efree);
    InfoPrint([0, 1], "outerSCF: Etot                        = ", scfDG.Etot);
    InfoPrint([0, 1], "outerSCF: norm(out-in)/norm(in) = ", scfDG.scfOuterNorm); 
    InfoPrint([0, 1], "outerSCF: Efree diff per atom   = ", scfDG.efreeDifPerAtom); 
    
    % Check for convergence
    if iter >= 2 ...
            && scfDG.scfOuterNorm < controlVar.scfOuterTolerance ...
            && scfDG.efreeDifPerAtom < controlVar.scfOuterEnergyTolerance
        % converged %
        InfoPrint([0, 1], ' Outer SCF is converged in %d steps ! \n', iter);
        isSCFConverged = true;
    end

    % Potential mixing for the outer SCF iteration. or no mixing at all anymore?
    % It seems that no mixing is the best.
        
    timeIterEnd = toc(timeIterStart);
    InfoPrint(0, ' Time for this SCF iteration = %8f [s]\n', timeIterEnd);
    
end  % -----------------  end of for iter -------------------------------


timeTotalEnd = toc(timeTotalStart);

InfoPrint([0, 1], '\nTotal time for all SCF iterations = %8f [s]\n', timeTotalEnd);

if isSCFConverged
    InfoPrint(0, ' Total number of outer SCF steps for SCF convergence = %d \n', iter-1);
else
    InfoPrint(0, ' Total number of outer SCF steps (SCF not converged) = %d \n', controlVar.scfOuterMaxIter);
end


% *********************************************************************
% Calculate the VDW contribution and the force
% *********************************************************************

timeForceStart = tic;
InfoPrint(0, '\n Computing forces using eigenvectors ...\n');
scfDG.hamDG = scfDG.hamDG.CalculateForce();
timeForceEnd = toc(timeForceStart);
InfoPrint(0, ' Time for computing the force is %8f [s] \n\n', timeForceEnd);


% Calculate the VDW energy
scfDG.Evdw = 0;
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

% output atomic structure information as file
if scfDG.userOption.isOutputAtomStruct
    domain = scfDG.domain;
    atomList = scfDG.hamDG.atomList;
    save(scfDG.dataFileIO.atomStructDG, 'domain', 'atomList');
end

% output data as files
% NOTE: data over each element will be stored in different files with the 
% same prefix from scfDG.dataFileIO

% output the density on the uniform fine grid
if scfDG.userOption.isOutputDensity
    uniformGridFine = scfDG.hamDG.grid.uniformGridElemFine;  % cell
    density = scfDG.hamDG.density;  % cell
    fileNamePrefix = scfDG.dataFileIO.densityDG;
    for elemIdx = 1 : numElemTotal
        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        densityElem = density{elemIdx};
        uniformGridFineElem = uniformGridFine{elemIdx};
        save(fileName, 'elemIdx', 'uniformGridFineElem', 'densityElem');
    end
end

% output potential on the global domain with uniform fine grid
if scfDG.userOption.isOutputPotential    
    uniformGridFine = scfDG.hamDG.grid.uniformGridElemFine;  % cell
    vtot = scfDG.hamDG.vtot;  % cell
    fileNamePrefix = scfDG.dataFileIO.potentialDG;
    for elemIdx = 1 : numElemTotal
        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        vtotElem = vtot{elemIdx};
        uniformGridFineElem = uniformGridFine{elemIdx};
        save(fileName, 'elemIdx', 'uniformGridFineElem', 'vtotElem');
    end    
end

% output wavefunction in the extended element
if scfDG.userOption.isOutputWfnExtElem
    fileNamePrefix = scfDG.dataFileIO.wfnExtElem;
    for elemIdx = 1 : numElemTotal
        eigSol = scfDG.vecEigSol{elemIdx};
        dm = eigSol.hamKS.domain;
        gridposExtElem = UniformMesh(dm);
        wfnExtElem = eigSol.psi.wavefun;
        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        save(fileName, 'elemIdx', 'gridposExtElem', 'wfnExtElem');
    end 
end

% output potential in the extended element
if scfDG.userOption.isOutputPotExtElem
    fileNamePrefix = scfDG.dataFileIO.potExtElem;
    for elemIdx = 1 : numElemTotal
        eigSol = scfDG.vecEigSol{elemIdx};
        dm = eigSol.hamKS.domain;
        gridposExtElemFine = UniformMeshFine(dm);
        vtotExtElem = eigSol.hamKS.vtot;
        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        save(fileName, 'elemIdx', 'gridposExtElemFine', 'vtotExtElem');
    end 
end

% output adaptive local basis (ALB) over uniform fine grid in each element
if scfDG.userOption.isOutputALBElemUniform
    hamDG = scfDG.hamDG;
    uniformGridFine = hamDG.grid.uniformGridElemFine;  % cell

    transferMat = hamDG.grid.LGLToUniformMatFine;    
    numGridOld = hamDG.grid.numLGLGridElem;
    numGridNew = hamDG.grid.numUniformGridElemFine;

    fileNamePrefix = scfDG.dataFileIO.albElemUniform;
    for elemIdx = 1 : numElemTotal
        basisLGL = hamDG.basisLGL{elemIdx};
        numBasis = size(basisLGL, 2);
        basisUniformFine = zeros(prod(numGridNew), numBasis);
        % interpolation from LGL grid to uniform fine grid
        for g = 1 : numBasis
            basisUniformFine(:, g) = InterpByTransferMat(transferMat, ...
                numGridOld, numGridNew, basisLGL(:, g));
        end
        uniformGridFineElem = uniformGridFine{elemIdx};

        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        save(fileName, 'elemIdx', 'uniformGridFineElem', 'basisUniformFine');
    end
end

% output adaptive local basis (ALB) over LGL grid in each element
if scfDG.userOption.isOutputALBElemLGL    
    hamDG = scfDG.hamDG;
    fileNamePrefix = scfDG.dataFileIO.albElemLGL;
    for elemIdx = 1 : numElemTotal
        gridposElemLGL = hamDG.grid.LGLGridElem{elemIdx};
        albElemLGL = hamDG.basisLGL{elemIdx}; 
        fileName = fileNamePrefix + "_" + num2str(elemIdx) + ".mat";
        save(fileName, 'elemIdx', 'gridposElemLGL', 'albElemLGL');
    end 
end


end