classdef SCFDG
    % SCFDG Self consistent iteration using the DG method    
    %
    %    scfDG = SCFDG() returns an empty SCFDG object.
    %
    %    scfDG = SCFDG(esdfParam, hamDG, vecEigSol, ptable) returns a 
    %    SCFDG object with respect to ESDFInputParam object esdfParam, 
    %    HamiltonianDG object hamDG, a cell vecEigSol containing 
    %    EigenSolverKS object over each extended element and a PeriodTable 
    %    object ptable.
    %
    %    See also HamiltonianDG, EigenSolverKS, PeriodTable, ESDFInputParam.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = private)
        domain
        hamDG
        vecEigSol 
        fft

        numElem
        numUnusedState
        XCType
        VDWType
        bufferSize
        
        controlVar        
        scfTotalInnerIter
        scfInnerNorm
        scfOuterNorm
        efreeDifPerAtom 
        
        Tbeta
        Tsigma
        vtotLGLSave
        EfreeHarris
        Efree
        Etot
        Ekin
        Ehart
        Ecor
        Exc
        Evdw
        EVxc
        Eself
        EIonSR
        Eext
        fermi
        forceVdw

        PeriodicUniformToLGLMat
        PeriodicUniformFineToLGLMat

        periodicPotential
        smearing
        mix                

        PWSolver
        CheFSIPW
        PPCGsbSize
        
        DGSolver
        
        userOption
        dataFileIO
    end
    
    methods (Access = public)
        function scfDG = SCFDG(varargin)
            scfDG.controlVar = struct(...
                'eigMinTolerance', [], ...
                'eigTolerance', [], ...
                'eigMinIter', [], ...
                'eigMaxIter', [], ...
                'scfInnerTolerance', [], ...
                'scfInnerMinIter', [], ...
                'scfInnerMaxIter', [], ...
                'scfOuterEnergyTolerance', [], ...
                'scfOuterTolerance', [], ...
                'scfOuterMinIter', [], ...
                'scfOuterMaxIter', [], ...
                'SVDBasisTolerance', [], ...
                ...
                'isPWeigToleranceDynamic', [] ...
                ...
                );
            
            scfDG.mix = struct(...
                'mixMaxDim', [], ...
                'mixType', [], ...
                'mixStepLength', [], ...
                'mixVariable', [], ...
                'mixOuterSave', [], ...
                'mixInnerSave', [], ...
                'dfOuterMat', [], ...
                'dvOuterMat', [], ...
                'dfInnerMat', [], ...
                'dvInnerMat', [] ...
                );
                        
            scfDG.periodicPotential = struct(...
                'distancePeriodize', [], ...
                'vBubble', [] ...
                );                
                                    
            scfDG.smearing = struct(...
                'SmearingScheme', "", ...
                'MPsmearingOrder', [] ...
                );
            
            scfDG.CheFSIPW = struct(...
                'firstFilterOrder', [], ...
                'firstCycleNum', [], ...
                'generalFilterOrder', [], ...
                'isApplyWfnEcut', [] ...
                );

            scfDG.userOption = struct(...
                'isOutputDensity', [], ...
                'isOutputPotenital', [], ...
                'isOutputAtomStruct', [], ...
                'isOutputALBElemUniform', [], ...
                'isOutputALBElemLGL', [], ...
                'isOutputWfnExtElem', [], ...
                'isOutputPotExtElem', [], ...
                ...
                'isPeriodizePotential', [], ...
                'isCalculateAPosterioriEachSCF', [] ...  % TODO
                );

            scfDG.dataFileIO = struct(...
                'densityDG', "", ...
                'potentialDG', "", ...
                'atomStructDG', "", ...
                'albElemUniform', "", ...
                'albElemLGL', "", ...
                'wfnExtElem', "", ...
                'potExtElem', "", ...
                ...
                'restartWfn', "", ...
                'restartDensity', "" ...
                );    
               
            switch (nargin)
                case 0 
                    scfDG.domain = Domain();
                    return;
                case 4
                    esdfParam = varargin{1};
                    hamDG = varargin{2};
                    vecEigSol = varargin{3};
                    ptable = varargin{4};
                    
                    fft = hamDG.fft;    
                    scfDG = Setup(scfDG, esdfParam, ...
                        hamDG, vecEigSol, fft, ptable);
                otherwise
                    error('Wrong number of arguments');
            end
        end
    end
    
end