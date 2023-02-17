function esdfParam = ESDFReadInput(inputFile)
% ESDFREADINPUT Input interface. Initialize ESDF parameters by inputFile
%    and default values, save data into ESDFInputParam object esdfParam.
%
%    See also ESDFInputParam, esdf_get, esdf_block.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


esdfParam =  ESDFInputParam();

% initialize and check parameters in inputFile
esdf_init(inputFile);


% ************************************************************************
%                            User Options
% ************************************************************************

% general options
esdfParam.userOption.general.isUseAtomDensity        = esdf_get( "Use_Atom_Density", 0 );
esdfParam.userOption.general.isUseVLocal             = esdf_get( "Use_VLocal", 0 );
esdfParam.userOption.general.isPWeigTolDynamic       = esdf_get( "Eig_Tolerance_Dynamic", 0 );
esdfParam.userOption.general.isRestartDensity        = esdf_get( "Restart_Density", 0 );
esdfParam.userOption.general.isRestartWfn            = esdf_get( "Restart_Wfn", 0 );

% PW
esdfParam.userOption.PW.isOutputWavefun              = esdf_get( "Output_Wavefun_PW", 0 );
esdfParam.userOption.PW.isOutputDensity              = esdf_get( "Output_Density_PW", 0 );
esdfParam.userOption.PW.isOutputPotential            = esdf_get( "Output_Potential_PW", 0 );
esdfParam.userOption.PW.isOutputAtomStruct           = esdf_get( "Output_Atom_Struct_PW", 0 );

% DG
esdfParam.userOption.DG.isOutputDensity               = esdf_get( "Output_Density_DG", 0 );
esdfParam.userOption.DG.isOutputPotential             = esdf_get( "Output_Potential_DG", 0 );
esdfParam.userOption.DG.isOutputAtomStruct            = esdf_get( "Output_Atom_Struct_DG", 0 );
esdfParam.userOption.DG.isOutputALBElemUniform        = esdf_get( "Output_ALB_Elem_Uniform", 0 );
esdfParam.userOption.DG.isOutputALBElemLGL            = esdf_get( "Output_ALB_Elem_LGL", 0 );
esdfParam.userOption.DG.isOutputWfnExtElem            = esdf_get( "Output_Wfn_ExtElem", 0 );
esdfParam.userOption.DG.isOutputPotExtElem            = esdf_get( "Output_Pot_ExtElem", 0 );
esdfParam.userOption.DG.isPeriodizePotential          = esdf_get( "Periodize_Potential", 0 );
esdfParam.userOption.DG.isCalculateAPosterioriEachSCF = esdf_get( "Calculate_APosteriori_Each_SCF", 0 );



% ************************************************************************
%                          Basic Parameters
% ************************************************************************

% -------------------------  Domain  -----------------------------------

[block_data, nlines] = esdf_block('Super_Cell');
if nlines ~= 0
    escell = sscanf(block_data(1), '%f');
    esdfParam.basic.domain.length(1) = escell(1);
    esdfParam.basic.domain.length(2) = escell(2);
    esdfParam.basic.domain.length(3) = escell(3);
else
    error('Super_Cell cannot be found.');
end
esdfParam.basic.domain.posStart = [0, 0, 0];

% NOTE: instead of assigning grid size, use ecut to determine the grid size 
% in global domain for wavefunction, density, and also other quantities in
% the local LGL domain. So domain.numGrid is not specified here.
% Grid of domain for PWDFT is combined into DG parameters part with 
% element size [1, 1, 1]. See the following code.

    
% --------------------------  Atoms  ---------------------------------

esdfParam.basic.atomList = Atom.empty();

esdfParam.basic.numAtomType = esdf_get('Atom_Types_Num', 0);
if esdfParam.basic.numAtomType == 0
    error('Atom_Types_Num cannot be found.');
end
type = esdf_get('Atom_Type', 0);  % type may be a vector
if type == 0
    error('Atom_Type cannot be found.');
end
esdfParam.basic.atomTypes = type;

% Read atomic coordinates
dm = esdfParam.basic.domain;
for i = 1 : esdfParam.basic.numAtomType
    
    % FIXME IMPORTANT. The "mass" parameter is removed from the reading
    % list. Mass can be obtained later with periodtable structure and the
    % atom type. 
    % NOTE that the mass in PeriodTable is in atomic mass unit (amu), but
    % the mass in atom vector is in atomic unit (au)
    
    [block_data, numAtom] = esdf_block('Atom_Bohr', i);
    % Cartesion coordinate (in the unit of Bohr)
    if numAtom ~= 0
        for j = 1 : numAtom
            pos = sscanf(block_data(j), '%f');
            pos = pos';  % turn column vector into row vector
            % force atom to be centered around (0, 0, 0)
            pos = pos - round( pos ./ dm.length ) .* dm.length;
            esdfParam.basic.atomList(end+1) = Atom(type(i), pos, [0, 0, 0], [0, 0, 0]);
        end
    else
    [block_data, numAtom] = esdf_block('Atom_Ang', i);
    % Cartesian coordinate (in the unit of angstrom)
    if numAtom ~= 0
        for k = 1 : numAtom
            pos = sscanf(block_data(k), '%f');
            pos = pos';  % turn column vector into row vector
            % unit conversion, from Ang to Bohr
            pos = pos ./ au2angDef();
            
            % force atom to be centered around (0, 0, 0)
            pos = pos - round( pos ./ dm.length ) .* dm.length;
            esdfParam.basic.atomList(end+1) = Atom(type(i), pos, [0, 0, 0], [0, 0, 0]);
        end
    else
    [block_data, numAtom] = esdf_block('Atom_Red', i);
    % Reduce coordinate (in the unit of Super_Cell)
    if numAtom ~= 0
        for l = 1 : numAtom
            pos = sscanf(block_data(l), '%f');
            pos = pos';  % turn column vector into row vector
            % unit conversion, from Reduced to Bohr
            pos = pos .* dm.length;
            
            % force atom to be centered around (0, 0, 0)
            pos = pos - round( pos ./ dm.length ) .* dm.length;
            esdfParam.basic.atomList(end+1) = Atom(type(i), pos, [0, 0, 0], [0, 0, 0]);
        end
    else
        error('Atomic coordinates cannot found for atom type: %d', type);
    end
    end
    end
end

        
% ------------------------  UPF_File  -----------------------------------

[block_data, nlines] = esdf_block('UPF_File');
if nlines ~= 0
    upfFile = block_data;
    upfFile = erase(upfFile, ' ');
    esdfParam.basic.upfFile = upfFile;
else
    esdfParam.basic.upfFile = [];
end

% -------------------------- Mixing -------------------------------------

% For metallic systems or small gapped semiconductors, mixStepLength is 
% often needed to be smaller than 0.1.  In such case, a better 
% preconditioner such as Kerker preconditioner can be helpful.

esdfParam.basic.mixMaxDim = esdf_get('Mixing_MaxDim', 9);

esdfParam.basic.mixType = esdf_get('Mixing_Type', 'anderson');
if (esdfParam.basic.mixType ~= "anderson" && ...
    esdfParam.basic.mixType ~= "kerker+anderson")
    error('Invalid mixing type.');
end

esdfParam.basic.mixVariable = esdf_get('Mixing_Variable', 'potential');
if (esdfParam.basic.mixVariable ~= "density" && ...
    esdfParam.basic.mixVariable ~= "potential")
    error('Invalid mixing variable');
end

esdfParam.basic.mixStepLength = esdf_get( "Mixing_StepLength", 0.8 );


% ------------------------  Others -------------------------------------

esdfParam.basic.ecutWavefunction = esdf_get( "Ecut_Wavefunction", 40.0 );
esdfParam.basic.densityGridFactor = esdf_get( "Density_Grid_Factor", 2.0 );

temperature = esdf_get("Temperature", 300.0);
esdfParam.basic.Tbeta = au2KDef() / temperature;


% Number of empty states for finite temperature calculation.
% This parameter must be larger than 0 for small gapped systems or
% relatively high temperature calculations.
esdfParam.basic.numExtraState  = esdf_get( "Extra_States",  0 );

% Some states for the planewave solver are unused in order to accelerate 
% the convergence rate of the eigensolver.
esdfParam.basic.numUnusedState = esdf_get( "Unused_States",  0 );

esdfParam.basic.extraElectron  = esdf_get( "Extra_Electron", 0 );

esdfParam.basic.pseudoType     = esdf_get("Pseudo_Type", "HGH"); 
esdfParam.basic.PWSolver       = esdf_get("PW_Solver", "LOBPCG"); 
esdfParam.basic.XCType         = esdf_get("XC_Type", "XC_LDA_XC_TETER93"); 
esdfParam.basic.VDWType        = esdf_get("VDW_Type", "None");
esdfParam.basic.smearingScheme = esdf_get("Smearing_Scheme", "FD");



% ************************************************************************
%                           Control Parameters
% ************************************************************************

esdfParam.control.scfInnerTolerance    = esdf_get( "SCF_Inner_Tolerance", 1e-4 );
esdfParam.control.scfInnerMinIter      = esdf_get( "SCF_Inner_MinIter",   1 );
esdfParam.control.scfInnerMaxIter      = esdf_get( "SCF_Inner_MaxIter",   1 );
esdfParam.control.scfOuterTolerance    = esdf_get( "SCF_Outer_Tolerance", 1e-6 );
esdfParam.control.scfOuterEnergyTolerance    = esdf_get( "SCF_Outer_Energy_Tolerance", 1e-4 );
esdfParam.control.scfOuterMinIter      = esdf_get( "SCF_Outer_MinIter",   3 );
esdfParam.control.scfOuterMaxIter      = esdf_get( "SCF_Outer_MaxIter",   30 );

% Default is no locking
esdfParam.control.eigTolerance         = esdf_get( "Eig_Tolerance", 1e-20 );
esdfParam.control.eigMinTolerance      = esdf_get( "Eig_Min_Tolerance", 1e-3 );
esdfParam.control.eigMinIter           = esdf_get( "Eig_MinIter",  2 );
esdfParam.control.eigMaxIter           = esdf_get( "Eig_MaxIter",  3 );
esdfParam.control.SVDBasisTolerance    = esdf_get( "SVD_Basis_Tolerance", 1e-6 );



% ************************************************************************
%                         PW parameters
% ************************************************************************

% PPCG
esdfParam.PW.PPCGsbSize = esdf_get( "PPCG_sbSize", 1 );

% Parameters related to Chebyshev Filtering in PWDFT
esdfParam.PW.CheFSI.firstFilterOrder   = esdf_get( "First_SCF_PWDFT_ChebyFilterOrder", 40 );
esdfParam.PW.CheFSI.firstCycleNum      = esdf_get( "First_SCF_PWDFT_ChebyCycleNum", 5 );
esdfParam.PW.CheFSI.generalFilterOrder = esdf_get( "General_SCF_PWDFT_ChebyFilterOrder", 35 );
esdfParam.PW.CheFSI.isApplyWfnEcut     = esdf_get( "PWDFT_Cheby_use_wfn_ecut_filt", 1 );



% ************************************************************************
%                           DG parameters
% ************************************************************************

% --------------------------  element  -----------------------------------

[block_data, nlines] = esdf_block("Element_Size");
if nlines ~= 0
    esdfParam.DG.numElem = sscanf(block_data, '%d', [1, 3]);
    esdfParam.isDGDFT = true;
else 
    esdfParam.isDGDFT = false;
    esdfParam.DG.numElem = [1, 1, 1];
end


if esdfParam.isDGDFT
    % Instead of grid size, use ecut to determine the number of grid points
    % in the local LGL domain.
    % The LGL grid factor does not need to be an integer
    esdfParam.DG.LGLGridFactor     = esdf_get( "LGL_Grid_Factor", 2.0 );

    esdfParam.DG.penaltyAlpha      = esdf_get( "Penalty_Alpha", 20.0 );
    
    % get the number of basis functions per elements
    % NOTE: ALB_Num_Element overwrites the parameter numALB later
    numALB = esdf_get("ALB_Num", 4);
    [block_data, nlines] = esdf_block("ALB_Num_Element");
    
    numElem = esdfParam.DG.numElem;
    numElemTotal = prod(numElem);
    if nlines ~= 0
        % use different number of ALB functions for each element
        esdfParam.DG.numALBElem = zeros(numElemTotal, 1);
        elemIdx = 0;
        for i = 1 : nlines
            data = sscanf(block_data(i), '%d');
            for j = 1 : length(data)
                elemIdx = elemIdx + 1;
                esdfParam.DG.numALBElem(elemIdx) = data(j);
            end
        end
        if elemIdx ~= numElemTotal
            error('The size of the number of ALB does not match the number of elements.');
        end            
    else
        % use the same number of ALB functions for each element.
        esdfParam.DG.numALBElem = numALB * ones(numElemTotal, 1);
    end

    % buffer size is number of elements in the buffer region along +(-) x(y,z) 
    % direction. Number of elements in extended element total will be 
    %    (2 * bufferSize + 1) ^ d
    % where d is number of dimensions that have been partitioned.    
    esdfParam.DG.bufferSize = esdf_get( "buffer_size", 1 );
    
    % Modification of the potential in the extended element    
    % Periodization of the external potential
    esdfParam.DG.distancePeriodize = [0, 0, 0];
    if esdfParam.userOption.DG.isPeriodizePotential
        [block_data, nlines] = esdf_block("Distance_Periodize");
        if nlines ~= 0
            distancePeriodize = sscanf(block_data, '%f');
            esdfParam.DG.distancePeriodize = reshape(distancePeriodize, 1, []);
        else
            % default value for DistancePeriodize
            for d = 1 : dimDef()
                if esdfParam.DG.numElem(d) == 1
                    esdfParam.DG.distancePeriodize(d) = 0;
                else
                    esdfParam.DG.distancePeriodize(d) = ... 
                        esdfParam.basic.domain.length(d) / numElem(d) * 0.5;
                end
            end
        end
    end  % end of modify the potential
    
    esdfParam.basic.DGSolver = esdf_get("DG_Solver", "eigs"); 

end

% --------------- Grid for PWDFT & DGDFT --------------------------------

% Choose the number of grid points
% NOTE: When applied to PWDFT, numElem = [1, 1, 1]
%
% The formula for the number of grid points along each dimension with 
% length L is 
%
%  1/2 K_max^2 = Ecut, with K_max = pi N_max / L;
%
% i.e., 
%
%  N_max = \frac{\sqrt{2 Ecut} * L}{pi}.
%
% The number of gird point along this dimension is chosen to be the largest
% even number bigger than N_max. The number of global grid points is also
% required to be divisible by the number of elements along that dimension.
%
% TODO Later the number of grid points can be improved to only contain the
% factor of 2, 3, and 5.

numElem = esdfParam.DG.numElem;
elemLength = esdfParam.basic.domain.length ./ numElem;
ecutWavefunction = esdfParam.basic.ecutWavefunction;
densityGridFactor = esdfParam.basic.densityGridFactor; 

% the number of grid is assumed to be at least an even number
numGridWavefunctionElem = ceil( ...
    sqrt(2 * ecutWavefunction) .* elemLength ./ pi ./ 2 ) .* 2;
numGridDensityElem = ceil( ...
    numGridWavefunctionElem .* densityGridFactor / 2 ) .* 2;

esdfParam.DG.numGridWavefunctionElem = numGridWavefunctionElem;
esdfParam.DG.numGridDensityElem = numGridDensityElem;

% coarse grid
esdfParam.basic.domain.numGrid = numGridWavefunctionElem .* numElem;
% fine grid
esdfParam.basic.domain.numGridFine = numGridDensityElem .* numElem;

if esdfParam.isDGDFT
    esdfParam.DG.numGridLGL = ceil( ...
        numGridWavefunctionElem * esdfParam.DG.LGLGridFactor ); 
end



% ************************************************************************
%                        IO Data File Names
% ************************************************************************

% PWDFT output
esdfParam.dataFileIO.wavefunPW    = esdf_get( "Wavefun_Output_File_PW", "WAVEFUN_PW.mat" );
esdfParam.dataFileIO.densityPW    = esdf_get( "Density_Output_File_PW", "DENSITY_PW.mat" );
esdfParam.dataFileIO.potentialPW  = esdf_get( "Potential_Output_File_PW", "POTENTIAL_PW.mat" );
esdfParam.dataFileIO.atomStructPW = esdf_get( "Atom_Struct_Output_File_PW", "ATOMSTRUCT_PW.mat" );

% DGDFT output
% NOTE: The following input/output options are prefix of data file names
% since IO datas are stored in element-wise format, data over different 
% elements will be saved in different files with file name format 
%   prefix + "_" + num2str(elemIdx) + ".mat" 

esdfParam.dataFileIO.densityDG      = esdf_get( "Density_Output_File_DG", "DENSITY_DG" );
esdfParam.dataFileIO.potentialDG    = esdf_get( "Potential_Output_File_DG", "POTENTIAL_DG" );
esdfParam.dataFileIO.atomStructDG   = esdf_get( "Atom_Struct_Output_File_DG", "ATOMSTRUCT_DG" );
esdfParam.dataFileIO.albElemUniform = esdf_get( "ALB_Elem_Uniform_Output_File", "ALB_UNIFORM" );
esdfParam.dataFileIO.albElemLGL     = esdf_get( "ALB_Elem_LGL_Output_File", "ALB_LGL" );
esdfParam.dataFileIO.wfnExtElem     = esdf_get( "Wavefun_ExtElem_Output_File", "WAVEFUN_EXTELEM" );
esdfParam.dataFileIO.potExtElem     = esdf_get( "Potential_ExtElem_Output_File", "POTENTIAL_EXTELEM" );

% restart data file
if esdfParam.isDGDFT
    esdfParam.dataFileIO.restartWfn = esdf_get( "Restart_Wfn_File", esdfParam.dataFileIO.wfnExtElem );
    esdfParam.dataFileIO.restartDensity = esdf_get( "Restart_Density_File", esdfParam.dataFileIO.densityDG );
else
    esdfParam.dataFileIO.restartWfn = esdf_get( "Restart_Wfn_File", esdfParam.dataFileIO.wavefunPW );
    esdfParam.dataFileIO.restartDensity = esdf_get( "Restart_Density_File", esdfParam.dataFileIO.densityPW );
end


% ***********************************************************************
% Some remaining consistency checks
% ***********************************************************************

if esdfParam.basic.pseudoType == "HGH" && ...
        esdfParam.userOption.general.isUseAtomDensity == true
    error('HGH type pseudopotential cannot use atom density as the initial guess');
end

        
end