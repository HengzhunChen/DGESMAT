classdef SCF
    % SCF class for self-consistent field iteration over the global domain 
    %    or extended element.
    %
    %    scf = SCF() returns an empty SCF object.
    %
    %    scf = SCF(eigSol) returns a SCF object with respect to 
    %    EigenSolverKS object eigSol.
    %
    %    See also EigenSolverKS, HamiltonianKS, PeriodTable.

    %  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
    %                     Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = private)
        controlVar
        PWSolver
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
        
        CheFSI           % Chebyshev filtering subspace iteration
        isChebyInIonDyn  % Do the usual Chebyshev filtering schedule 
                         % or work in ion dynamic mode
        mixing
        restart
    end
    
    methods
        function scf = SCF(varargin)
            scf.controlVar = struct(...
                'eigMinTolerance', [], ...
                'eigTolerance', [], ...
                'eigMaxIter', [], ...
                'scfTolerance', [], ...
                'scfMaxIter', [], ...
                'scfPhiMaxIter', [], ...
                'scfPhiTolerance', [], ...
                'numUnusedState', [], ...
                'isEigToleranceDynamic', [], ...
                'isCalculateGradRho', [] ...  % need for GGA, meta-GGA and hybrid functional
                );
            
            scf.mixing = struct(...
                'mixMaxDim', [], ...
                'mixType', [], ...
                'mixStepLength', [], ...
                'dfMat', [], ...
                'dvMat', [], ...
                'mixSave', [] ...  % work array for the mixing variable in the inner iteration
                );
            
            scf.restart = struct(...
                'DensityFileName', "", ...
                'WfnFileName', "" ...
                );
            
            scf.CheFSI = struct(...
                'firstFilterOrder', [], ...
                'firstCycleNum', [], ...
                'generalFilterOrder', [], ...
                'isApplyWfnEcut', [] ...
                );

        
            switch (nargin)
                case 0            
                    scf.eigSol = EigenSolverKS();
                case 1
                    eigSol = varargin{1};
                    scf = Setup(scf, eigSol);
            end
        end              
    end
end        