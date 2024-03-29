classdef ESDFInputParam
    % Main class containing input parameters for the electronic structure 
    % calculation. ESDF means Electronic Structure Data Format.
    %
    %    esdfParam = ESDFInputParam() returns an ESDFInputParam object
    %    without initialization, which is done by ESDFReadInput().
    %
    %    This ESDF structure is designed to simplify and enhance the input 
    %    of data into electronic structure codes. 
    %    An important feature is the requirement that most inputs require 
    %    default settings to be supplied within the main program calling
    %    esdf. This means that rarely used variables will not clutter 
    %    everyday input files, and, even more usefully, "intelligence" 
    %    may be built into the main code as the defaults may be dependent 
    %    of other set variables. 
    %    Another important feature is the ability to define "physical" 
    %    values. This means that the input files need not depend on the 
    %    internal physical units used by the main program.
    %
    %    To add new input variable into esdf structure, 
    %    1. update the new key word in esdf_key()
    %    2. use esdf_get() or esdf_block() to read data from input file or
    %    set default value in ESDFReadInput(), and data will be passed to
    %    global variable esdfParam.
    %
    %    For more explanation of each variables, see also the comment and
    %    note in ESDFReadInput.
    %
    %    See also ESDFReadInput, ESDFPrintInput.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties ( SetAccess = public )
        basic
        control
        dataFileIO
        
        PW
        DG
                
        userOption
        devOption
        
        isDGDFT
        % This is NOT an input parameter, but records whether
        % DGDFT is performed
    end
  
    methods
        function esdfParam = ESDFInputParam()
            
            esdfParam.basic = struct(...
                'domain', Domain(), ...
                'atomList', Atom.empty(), ...
                'numAtomType', [], ...
                'atomTypes', [], ...
                'ecutWavefunction', [], ...
                'densityGridFactor', [], ...
                ...
                'Tbeta', [], ...
                'numExtraState', [], ...
                'numUnusedState', [], ...
                'extraElectron', [], ...
                ...
                'mixType', "", ...
                'mixVariable', "", ...
                'mixMaxDim', [], ...
                'mixStepLength', [], ...
                ...
                'pseudoType', "", ...
                'upfFile', [], ...
                'XCType', "", ...
                'VDWType', "", ...
                'smearingScheme', "", ...
                ...
                'PWSolver', "", ...
                'DGSolver', "" ...
                 );

             esdfParam.control = struct(...
                'eigMinTolerance', [],...
                'eigTolerance', [], ...
                'eigMinIter', [], ...
                'eigMaxIter', [], ...
                ...
                'scfInnerMinIter', [], ...
                'scfInnerMaxIter', [], ...
                'scfOuterMinIter', [], ...
                'scfOuterMaxIter', [], ...
                'scfInnerTolerance', [], ...
                'scfOuterTolerance', [], ...
                'scfOuterEnergyTolerance', [], ...
                ...
                'SVDBasisTolerance', [] ...
                );

            esdfParam.PW = struct(...
                'PPCGsbSize', [], ...
                'CheFSI', [] ...
                );
            
            esdfParam.PW.CheFSI = struct( ...
                'firstFilterOrder', [], ...
                'firstCycleNum', [], ...
                'generalFilterOrder', [], ...
                'isApplyWfnEcut', [] ...
                );

            esdfParam.DG = struct(...
                'numElem', [], ...
                'numGridWavefunctionElem', [], ...
                'numGridDensityElem', [], ...
                'numGridLGL', [], ...
                'numALBElem', [], ...
                'LGLGridFactor', [], ...
                'bufferSize', [], ...
                'penaltyAlpha', [], ...
                ...
                'distancePeriodize', [] ...
                );
                        
            esdfParam.userOption = struct(...
                'general', [], ...
                'PW', [], ...
                'DG', [] ...
                );

            esdfParam.userOption.general = struct(...
                'isRestartWfn', [], ...
                'isRestartDensity', [], ...
                ...
                'isUseAtomDensity', [], ...
                'isPWeigTolDynamic', [] ...
                );

            esdfParam.userOption.PW = struct(...
                'isOutputWavefun', [], ...
                'isOutputDensity', [], ...
                'isOutputPotential', [], ...
                'isOutputAtomStruct', [] ...
            );

            esdfParam.userOption.DG = struct(...
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

            esdfParam.dataFileIO = struct(...
                ...  % output file names for PWDFT           
                'wavefunPW', "", ...
                'densityPW', "", ...
                'potentialPW', "", ...
                'atomStructPW', "", ...
                ...  % output file names for DGDFT
                'densityDG', "", ...
                'potentialDG', "", ...
                'atomStructDG', "", ...
                'albElemUniform', "", ...
                'albElemLGL', "", ...
                'wfnExtElem', "", ...
                'potExtElem', "", ...
                ...  % restart file names for PWDFT or DGDFT
                'restartWfn', "", ...
                'restartDensity', "" ...
                );

            % Advanced options for developers
            esdfParam.devOption = struct(...
                ...
                );
            
        end
    end
end