classdef SCF
    % SCF class for self-consistent field iteration over the global domain 
    %    or extended element.
    %
    %    scf = SCF(esdfParam) returns an empty SCF object according to 
    %    ESDFInputParam object esdfParam.
    %
    %    scf = SCF(esdfParam, eigSol) returns a SCF object with respect to 
    %    ESDFInputParam object esdfParam, EigenSolverKS object eigSol.
    %
    %    See also EigenSolverKS, HamiltonianKS, PeriodTable, ESDFInputParam.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = private)
        controlVar
        eigSol
        
        % Physical parameters %
        Tbeta                    % Inverse of temperature in atomic unit
        Efree                    % Helmholtz free energy
        EfreeHarris              % Free energy from Harris functional
        Etot                     % Total energy
        Ekin                     % Kinetic energy
        Ehart                    % Hartree (correction) energy
        Ecor                     % Nonlinear correction energy
        Exc                      % Exchange-correlation energy
        EVdw                     % Van der Waals energy
        EVxc                     % Exchange-correlation potential energy
        Eself                    % Self energy due to the pseudopotential
        EIonSR                   % Short range repulsion energy for Gaussian charge
        Eext                     % External energy
        fermi                    % Fermi energy
        Efock                    % Hartree-Fock energy
        totalCharge              % Total number of computed electron charge
        
        % SCF variables %
        vtotNew
        scfNorm          % || V_{new}-V_{old} || / ||V_{old}||
        efreeDifPerAtom  % difference bwtween free energy and Harris free energy per atom        

        PWSolver
        PPCGsbSize       % parameter of PPCG, size of subblocks in PPCG
        CheFSI           % Chebyshev filtering subspace iteration

        mixing

        userOption
        dataFileIO
    end
    
    methods
        function scf = SCF(varargin)
            scf.controlVar = struct(...
                'eigMinTolerance', [], ...
                'eigTolerance', [], ...
                'eigMaxIter', [], ...
                'scfTolerance', [], ...
                'scfMaxIter', [], ...
                'isEigToleranceDynamic', [], ...
                'isCalculateGradRho', [] ...  % need for GGA, meta-GGA and hybrid functional
                );
            
            scf.mixing = struct(...
                'mixMaxDim', [], ...
                'mixType', [], ...
                'mixStepLength', [], ...
                'mixVariable', [], ...
                'dfMat', [], ...
                'dvMat', [], ...
                'mixSave', [] ...  % work array for the mixing variable in the inner iteration
                );
                        
            scf.CheFSI = struct(...
                'firstFilterOrder', [], ...
                'firstCycleNum', [], ...
                'generalFilterOrder', [], ...
                'isApplyWfnEcut', [] ...
                );

            scf.userOption = struct(...
                'isOutputWavefun', [], ...
                'isOutputDensity', [], ...
                'isOutputPotential', [], ...
                'isOutputAtomStruct', [] ...
                ...
                );

            scf.dataFileIO = struct(...
                'wavefunPW', "", ...
                'densityPW', "", ...
                'potentialPW', "", ...
                'atomStructPW', "", ...
                ...
                'restartWfn', "", ...
                'restartDensity', "" ...
                );

        
            switch (nargin)
                case 1
                    esdfParam = varargin{1};
                    eigSol = EigenSolverKS();
                    scf = Setup(scf, esdfParam, eigSol);
                case 2
                    esdfParam = varargin{1};
                    eigSol = varargin{2};
                    scf = Setup(scf, esdfParam, eigSol);
            end
        end              
    end
end        