function ESDFReadInput(inputFile)
% ESDFREADINPUT Input interface. Initialize ESDF parameters by inputFile
%               and default values, save data into global variable
%               esdfParam.
%
%    See also ESDFInputParam, esdf_get, esdf_block.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global esdfParam;

% initialize and check parameters in inputFile
esdf_init(inputFile);


% ************************************************************************
%                            User Options
% ************************************************************************

% general options
esdfParam.userOption.general.isUseAtomDensity        = esdf_get( "Use_Atom_Density", 0 );
esdfParam.userOption.general.isUseVLocal             = esdf_get( "Use_VLocal", 0 );
esdfParam.userOption.general.isRestartDensity        = esdf_get( "Restart_Density", 0 );
esdfParam.userOption.general.isRestartWfn            = esdf_get( "Restart_Wfn", 0 );
esdfParam.userOption.general.isOutputDensity         = esdf_get( "Output_Density", 0 );
esdfParam.userOption.general.isOutputStructInfo      = esdf_get( "Output_Struct_Info", 0 );
esdfParam.userOption.general.isPWeigTolDynamic       = esdf_get( "Eig_Tolerance_Dynamic", 0 );

% PW
esdfParam.userOption.PW.isOutputWfn                  = esdf_get( "Output_Wfn", 0 );

% DG
esdfParam.userOption.DG.isOutputWfnExtElem            = esdf_get( "Output_Wfn_ExtElem", 0 );
esdfParam.userOption.DG.isPotentialBarrier            = esdf_get( "Potential_Barrier",  0 );
esdfParam.userOption.DG.isPeriodizePotential          = esdf_get( "Periodize_Potential", 0 );
esdfParam.userOption.DG.isCalculateAPosterioriEachSCF = esdf_get( "Calculate_APosteriori_Each_SCF", 0 );
esdfParam.userOption.DG.isCalculateForceEachSCF       = esdf_get( "Calculate_Force_Each_SCF", 0 );

% ionDyn
esdfParam.userOption.ionDyn.isRestartPosition     = esdf_get( "Restart_Position", 0 );
esdfParam.userOption.ionDyn.isRestartVelocity     = esdf_get( "Restart_Velocity", 0 );
esdfParam.userOption.ionDyn.isOutputPosition      = esdf_get( "Output_Position", 1 );
esdfParam.userOption.ionDyn.isOutputVelocity      = esdf_get( "Output_Velocity", 1 );
esdfParam.userOption.ionDyn.isOutputXYZ           = esdf_get( "Output_XYZ", 1 );



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

% Read position from lastPos.out into esdfParam.atomList[i].pos
if esdfParam.userOption.ionDyn.isRestartPosition
    InfoPrint([0, 1], '\nRead in atomic position from lastPos.out, \n');
    InfoPrint([0, 1], 'override the atomic positions read from the input file.\n');
    % read atom position from lastPos.out
    inFid = fopen('lastPos.out', 'r');
    atomposRead = fscanf(inFid, '%f');
    for i = 1 : length(esdfParam.basic.atomList)-1
        esdfParam.basic.atomList(i).pos = atomposRead(3*i: 3*i+3);
    end
    fclose(inFid);
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

% For metallic systems or small gapped semiconductors, mixStepLength is 
% often needed to be smaller than 0.1.  In such case, a better 
% preconditioner such as Kerker preconditioner can be helpful.

% ------------------------  Others -------------------------------------

esdfParam.basic.ecutWavefunction     = esdf_get( "Ecut_Wavefunction", 40.0 );
esdfParam.basic.densityGridFactor    = esdf_get( "Density_Grid_Factor", 2.0 );

temperature = esdf_get("Temperature", 300.0);
esdfParam.basic.Tbeta = au2KDef() / temperature;

esdfParam.basic.numExtraState  = esdf_get( "Extra_States",  0 );
esdfParam.basic.numUnusedState = esdf_get( "Unused_States",  0 );
esdfParam.basic.extraElectron  = esdf_get( "Extra_Electron", 0 );


esdfParam.basic.pseudoType      = esdf_get("Pseudo_Type", "HGH"); 
esdfParam.basic.PWSolver        = esdf_get("PW_Solver", "LOBPCG"); 
esdfParam.basic.XCType          = esdf_get("XC_Type", "XC_LDA_XC_TETER93"); 
esdfParam.basic.VDWType         = esdf_get("VDW_Type", "None");
esdfParam.basic.smearingScheme  = esdf_get("Smearing_Scheme", "FD");



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
esdfParam.PW.CheFSI.firstFilterOrder   = esdf_get("First_SCF_PWDFT_ChebyFilterOrder", 40 );
esdfParam.PW.CheFSI.firstCycleNum      = esdf_get("First_SCF_PWDFT_ChebyCycleNum", 5 );
esdfParam.PW.CheFSI.generalFilterOrder = esdf_get("General_SCF_PWDFT_ChebyFilterOrder", 35 );
esdfParam.PW.CheFSI.isApplyWfnEcut     = esdf_get("PWDFT_Cheby_use_wfn_ecut_filt", 1 );



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
    % The LGL grid facotr does not need to be an integer
    esdfParam.DG.LGLGridFactor     = esdf_get( "LGL_Grid_Factor", 2.0 );
    esdfParam.DG.GaussInterpFactor = esdf_get( "Gauss_Interp_Factor", 4.0 );
    esdfParam.DG.GaussSigma        = esdf_get( "Gauss_Sigma", 0.001 );
    esdfParam.DG.penaltyAlpha      = esdf_get( "Penalty_Alpha", 20.0 );
    
    % get the number of basis functions per elements
    % NOTE: ALB_Num_Element overwrites the parameter numALB later
    numALB = esdf_get("ALB_Num", 4);
    [block_data, sizeALBElem] = esdf_block("ALB_Num_Element");
    
    numElem = esdfParam.DG.numElem;
    if sizeALBElem ~= 0
        % use different number of ALB functions for each element
        if sizeALBElem ~= prod(numElem)
            error('The size of the number of ALB does not match the number of elements.');
        end
        for k = 1 : numElem(3)
            for j = 1 : numElem(2)
                for i = 1 : numElem(1)
                    esdfParam.DG.numALBElem(i, j, k) = sscanf( block_data( ... 
                    i + (j-1)*numElem(1) + (k-1)*numElem(1)*numElem(2), '%d' ));
                end
            end
        end
    else
        % use the same number of ALB functions for each element.
        esdfParam.DG.numALBElem = numALB * ones(numElem(1), numElem(2), numElem(3));
    end
    
    % Modification of the potential in the extended element
    % FIXME the potential barrier is now obsolete.
    esdfParam.DG.potentialBarrierW    = esdf_get( "Potential_Barrier_W", 2.0 );
    esdfParam.DG.potentialBarrierS    = esdf_get( "Potential_Barrier_S", 0.0 );
    esdfParam.DG.potentialBarrierR    = esdf_get( "Potential_Barrier_R", 5.0 );
    
    % Periodization of the external potential
    esdfParam.DG.distancePeriodize = [0, 0, 0];
    if esdfParam.userOption.DG.isPeriodizePotential
        [block_data, nlines] = esdf_block("Distance_Periodize");
        if nlines ~= 0
            esdfParam.DG.distancePeriodize = sscanf(block_data, '%f');
        else
            % default value for DistancePeriodize
            for d = 1 : dimDef()
                if esdfParam.DG.numElem(d) == 1
                    esdfParam.DG.distancePeriodize(d) = 0;
                else
                    esdfParam.DG.distancePeriodize(d) = ... 
                        esdfParam.basic.domain.length(d) / esdfParam.DG.numElem(d) * 0.5;
                end
            end
        end
    end
    % end of modify the potential
    
    
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

elemLength = esdfParam.basic.domain.length ./ esdfParam.DG.numElem;

% the number of grid is assumed to be at least an even number
esdfParam.DG.numGridWavefunctionElem = ceil( ...
    sqrt(2 * esdfParam.basic.ecutWavefunction) .* elemLength ./ pi ./ 2 ) .* 2;
esdfParam.DG.numGridDensityElem = ceil( ...
    esdfParam.DG.numGridWavefunctionElem .* esdfParam.basic.densityGridFactor / 2 ) .* 2;

% coarse grid
esdfParam.basic.domain.numGrid = esdfParam.DG.numGridWavefunctionElem .* esdfParam.DG.numElem;
% fine grid
esdfParam.basic.domain.numGridFine = esdfParam.DG.numGridDensityElem .* esdfParam.DG.numElem;

if esdfParam.isDGDFT
    esdfParam.DG.numGridLGL = ceil( esdfParam.DG.numGridWavefunctionElem * esdfParam.DG.LGLGridFactor ); 
end


% ----------------------- CheFSI DG parameters ---------------------------

if esdfParam.isDGDFT

% First SCF step parameters
esdfParam.DG.CheFSI.firstFilterOrder = esdf_get( "First_SCFDG_ChebyFilterOrder", 60 );
esdfParam.DG.CheFSI.firstCycleNum    = esdf_get( "First_SCFDG_ChebyCycleNum", 5 );

% Second stage parameters
esdfParam.DG.CheFSI.secondOuterIter   = esdf_get( "Second_SCFDG_ChebyOuterIter", 3 );
esdfParam.DG.CheFSI.secondFilterOrder = esdf_get( "Second_SCFDG_ChebyFilterOrder", 60 );
esdfParam.DG.CheFSI.secondCycleNum    = esdf_get( "Second_SCFDG_ChebyCycleNum", 3);

% General SCF step parameters
esdfParam.DG.CheFSI.generalFilterOrder = esdf_get( "General_SCFDG_ChebyFilterOrder", 60);
esdfParam.DG.CheFSI.generalCycleNum    = esdf_get( "General_SCFDG_ChebyCycleNum", 1);

% -------------- Complementary subspace strategy in CheFSI --------------

esdfParam.DG.CheFSI.isUseCompSubspace = esdf_get("SCFDG_use_CheFSI_complementary_subspace", 0);

esdfParam.DG.CheFSI.CompSubspace.nStates = ...
    esdf_get("SCFDG_complementary_subspace_nstates", fix(esdfParam.basic.numExtraState/20.0 + 0.5) );
esdfParam.DG.CheFSI.CompSubspace.ionIterRegularChebyFreq = ...
    esdf_get("SCFDG_CS_ioniter_regular_Cheby_freq", 20 );

esdfParam.DG.CheFSI.CompSubspace.biggerGridDimFac = ...
    esdf_get("SCFDG_CS_bigger_grid_dim_fac", 1 );

% Inner LOBPCG related options
esdfParam.DG.CheFSI.CompSubspace.lobpcgIter = ...
    esdf_get("SCFDG_complementary_subspace_inner_LOBPCGiter", 15);
esdfParam.DG.CheFSI.CompSubspace.lobpcgTol = ...
    esdf_get("SCFDG_complementary_subspace_inner_LOBPCGtol", 1e-8);

% Inner CheFSI related options
esdfParam.DG.CheFSI.CompSubspace.isHMatTopStatesUseCheby = ...
    esdf_get("SCFDG_complementary_subspace_use_inner_Cheby", 1);
esdfParam.DG.CheFSI.CompSubspace.HMatFilterOrder =  ...
    esdf_get("SCFDG_complementary_subspace_inner_Chebyfilterorder", 5);
esdfParam.DG.CheFSI.CompSubspace.HMatCycleNum = ...
    esdf_get("SCFDG_complementary_subspace_inner_Chebycyclenum", 3);

end



% ************************************************************************
%                        Parameters for Hybird
% ************************************************************************

esdfParam.hybrid.scfPhiMaxIter   = esdf_get( "SCF_Phi_MaxIter",     10 );
esdfParam.hybrid.scfPhiTolerance = esdf_get( "SCF_Phi_Tolerance",   1e-6 );

esdfParam.hybrid.MixType = esdf_get("Hybrid_Mixing_Type", "nested");
if esdfParam.hybrid.MixType ~= "nested" && ...
   esdfParam.hybrid.MixType ~= "scdiis" && ...
   esdfParam.hybrid.MixType ~= "pcdiis"
    error('Invalid hybrid mixing type.');
end

esdfParam.hybrid.isHybridACETwicePCDIIS = esdf_get( "Hybrid_ACE_Twice_PCDIIS", 1 );
esdfParam.hybrid.isHybridACE            = esdf_get( "Hybrid_ACE", 1 );
esdfParam.hybrid.isHybridActiveInit     = esdf_get( "Hybrid_Active_Init", 0 );
esdfParam.hybrid.isHybridDF             = esdf_get( "Hybrid_DF", 0 );

esdfParam.hybrid.DFType = esdf_get( "Hybrid_DF_Type", "QRCP" );
if esdfParam.hybrid.DFType ~= "QRCP" && ...
   esdfParam.hybrid.DFType ~= "Kmeans" && ...
   esdfParam.hybrid.DFType ~= "Kmeans+QRCP"
    error('Invalid ISDF type.');
end

esdfParam.hybrid.DFKmeansWFType = esdf_get( "Hybrid_DF_Kmeans_WF_Type", "Add" );
if esdfParam.hybrid.DFKmeansWFType ~= "Add" && ...
   esdfParam.hybrid.DFKmeansWFType ~= "Multi"
    error('Invalid Kmeans WF type.');
end

esdfParam.hybrid.DFKmeansWFAlpha     = esdf_get( "Hybrid_DF_Kmeans_WF_Alpha", 2.0 ); % 0.5 1.0 1.5 2.0 2.5 3.0 4.0
esdfParam.hybrid.DFKmeansTolerance   = esdf_get( "Hybrid_DF_Kmeans_Tolerance", 1e-3 );
esdfParam.hybrid.DFKmeansMaxIter     = esdf_get( "Hybrid_DF_Kmeans_MaxIter", 99 );
esdfParam.hybrid.DFNumMu             = esdf_get( "Hybrid_DF_Num_Mu", 6.0 );
esdfParam.hybrid.DFNumGaussianRandom = esdf_get( "Hybrid_DF_Num_GaussianRandom", 2.0 );
esdfParam.hybrid.DFTolerance         = esdf_get( "Hybrid_DF_Tolerance", 1e-20 );

esdfParam.hybrid.MDscfPhiMaxIter     = esdf_get( "MD_SCF_Phi_MaxIter", esdfParam.hybrid.scfPhiMaxIter  );
% This is used in DGDFT for energy based SCF
esdfParam.hybrid.MDscfOuterMaxIter   = esdf_get( "MD_SCF_Outer_MaxIter",  esdfParam.control.scfOuterMaxIter ); 

esdfParam.hybrid.exxDivergenceType   = esdf_get( "EXX_Divergence_Type", 1 );



% ************************************************************************
%                    Parameters for Ionic Motion
% ************************************************************************

% Both for geometry optimization and molecular dynamics
% The default is 0, which means that only static claculation
esdfParam.ionDyn.ionMaxIter = esdf_get("Ion_Max_Iter", 0);
esdfParam.ionDyn.ionMove    = esdf_get("Ion_Move", "");

% geometry optimization
esdfParam.ionDyn.geoOptMaxForce = esdf_get("Geo_Opt_Max_Force", 0.001);

% NLCG related parameters
esdfParam.ionDyn.geoOptNLCGsigma = esdf_get("Geo_Opt_NLCG_Sigma", 0.02);

% FIRE related parameters
esdfParam.ionDyn.fireNmin = esdf_get("FIRE_Nmin", 5);  % Compare with LAMMPS
esdfParam.ionDyn.firedt = esdf_get("FIRE_Time_Step", 41.3413745758);  % usually between 0.1-1fs 
esdfParam.ionDyn.fireAtomicMass = esdf_get("FIRE_Atomic_Mass", 4.0);  % Compare with LAMMPS

% Molecular dynamics
esdfParam.ionDyn.ionTemperature        = esdf_get("Ion_Temperature", 300.0);
esdfParam.ionDyn.TbetaIonTemperature   = au2KDef() / esdfParam.ionDyn.ionTemperature;

esdfParam.ionDyn.MDTimeStep              = esdf_get("MD_Time_Step", 40.0);
esdfParam.ionDyn.MDExtrapolationType     = esdf_get("MD_Extrapolation_Type", "linear"); 
esdfParam.ionDyn.MDExtrapolationVariable = esdf_get("MD_Extrapolation_Variable", "density");
esdfParam.ionDyn.qMass                   = esdf_get("Thermostat_Mass", 85000.0);
esdfParam.ionDyn.langevinDamping         = esdf_get("Langevin_Damping", 0.01);
esdfParam.ionDyn.kappaXLBOMD             = esdf_get("kappa_XLBOMD", 1.70);

% Energy based SCF convergence for MD: currently used in DGDFT only
esdfParam.ionDyn.MDscfEnergyCriteriaEngageIonIter = ...
    esdf_get("MD_SCF_energy_criteria_engage_ioniter", esdfParam.ionDyn.ionMaxIter + 1); 
esdfParam.ionDyn.MDscfEtotDiffTol = ...
    esdf_get("MD_SCF_Etot_diff", esdfParam.control.scfOuterEnergyTolerance);
esdfParam.ionDyn.MDscfEbandDiffTol = ...
    esdf_get("MD_SCF_Eband_diff", esdfParam.control.scfOuterEnergyTolerance);



% ***********************************************************************
% Some remaining consistency checks
% ***********************************************************************

if esdfParam.basic.pseudoType == "HGH" && ...
        esdfParam.userOption.general.isUseAtomDensity == true
    error('HGH type pseudopotential cannot use atom density as the initial guess');
end

        
end