function esdf_key()
% ESDF_KEY Module to hold keyword list, used to initialize global variable
%    kw_label.
%
%    This must be updated as new keywords are brought into existence.
%
%    See also esdf_init, ESDFReadInput, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


global kw_label;

% ************************************************************************
%                            User Options
% ************************************************************************

i = 1;
kw_label(i) = "use_atom_density";

i = i + 1;
kw_label(i) = "use_vlocal";

i = i + 1;
kw_label(i) = "eig_tolerance_dynamic";


i = i + 1;
kw_label(i) = "restart_density";

i = i + 1;
kw_label(i) = "restart_wfn";


i = i + 1;
kw_label(i) = "output_wavefun_pw";

i = i + 1;
kw_label(i) = "output_density_pw";

i = i + 1;
kw_label(i) = "output_potential_pw";

i = i + 1;
kw_label(i) = "output_atom_struct_pw";


i = i + 1;
kw_label(i) = "output_density_dg";

i = i + 1;
kw_label(i) = "output_potential_dg";

i = i + 1;
kw_label(i) = "output_atom_struct_dg";

i = i + 1;
kw_label(i) = "output_alb_elem_lgl";

i = i + 1;
kw_label(i) = "output_alb_elem_uniform";

i = i + 1;
kw_label(i) = "output_wfn_extelem";

i = i + 1;
kw_label(i) = "output_pot_extelem";


i = i + 1;
kw_label(i) = "periodize_potential";

i = i + 1;
kw_label(i) = "calculate_aposteriori_each_scf";



% ************************************************************************
%                        IO Data File Names 
% ************************************************************************

i = i + 1;
kw_label(i) = "wavefun_output_file_pw";

i = i + 1;
kw_label(i) = "density_output_file_pw";

i = i + 1;
kw_label(i) = "potential_output_file_pw";

i = i + 1;
kw_label(i) = "atom_struct_output_file_pw";


i = i + 1;
kw_label(i) = "density_output_file_dg";

i = i + 1;
kw_label(i) = "potential_output_file_dg";

i = i + 1;
kw_label(i) = "atom_struct_output_file_dg";

i = i + 1;
kw_label(i) = "alb_elem_uniform_output_file";

i = i + 1;
kw_label(i) = "alb_elem_lgl_output_file";

i = i + 1;
kw_label(i) = "wavefun_extelem_output_file";

i = i + 1;
kw_label(i) = "potential_extelem_output_file";


i = i + 1;
kw_label(i) = "restart_wfn_file";

i = i + 1;
kw_label(i) = "restart_density_file";



% ************************************************************************
%                          Basic Parameters
% ************************************************************************

i = i + 1;
kw_label(i) = "atom_types_num";

i = i + 1;
kw_label(i) = "atom_type";

i = i + 1;
kw_label(i) = "atom_coord";

i = i + 1;
kw_label(i) = "atom_bohr";

i = i + 1;
kw_label(i) = "atom_ang";

i = i + 1;
kw_label(i) = "atom_red";

i = i + 1;
kw_label(i) = "super_cell";

i = i + 1;
kw_label(i) = "upf_file";

i = i + 1;
kw_label(i) = "mixing_type";

i = i + 1;
kw_label(i) = "mixing_variable";

i = i + 1;
kw_label(i) = "mixing_maxdim";

i = i + 1;
kw_label(i) = "mixing_steplength";

i = i + 1;
kw_label(i) = "extra_states";

i = i + 1;
kw_label(i) = "pseudo_type";

i = i + 1;
kw_label(i) = "pw_solver";

i = i + 1;
kw_label(i) = "xc_type";

i = i + 1;
kw_label(i) = "vdw_type";

i = i + 1;
kw_label(i) = "smearing_scheme";

i = i + 1;
kw_label(i) = "temperature";

i = i + 1;
kw_label(i) = "extra_electron";

i = i + 1;
kw_label(i) = "unused_states";



% ************************************************************************
%                           Control Parameters
% ************************************************************************

i = i + 1;
kw_label(i) = "eig_tolerance";

i = i + 1;
kw_label(i) = "eig_min_tolerance";

i = i + 1;
kw_label(i) = "eig_miniter";

i = i + 1;
kw_label(i) = "eig_maxiter";

i = i + 1;
kw_label(i) = "scf_inner_tolerance";

i = i + 1;
kw_label(i) = "scf_outer_tolerance";

i = i + 1;
kw_label(i) = "scf_outer_energy_tolerance";

i = i + 1;
kw_label(i) = "svd_basis_tolerance";

i = i + 1;
kw_label(i) = "scf_inner_miniter";

i = i + 1;
kw_label(i) = "scf_inner_maxiter";

i = i + 1;
kw_label(i) = "scf_outer_miniter";

i = i + 1;
kw_label(i) = "scf_outer_maxiter";



% ************************************************************************
%                         PW parameters
% ************************************************************************

i = i + 1;
kw_label(i) = "ppcg_sbsize";

% Inputs related to Chebyshev polynomial Filtered SCF iterations for PWDFT    
i = i + 1;
kw_label(i) = "first_scf_pwdft_chebyfilterorder";

i = i + 1;
kw_label(i) = "first_scf_pwdft_chebycyclenum";

i = i + 1;
kw_label(i) = "general_scf_pwdft_chebyfilterorder";

i = i + 1;
kw_label(i) = "pwdft_cheby_use_wfn_ecut_filt";



% ************************************************************************
%                           DG parameters
% ************************************************************************


i = i + 1;
kw_label(i) = "alb_num";

i = i + 1;
kw_label(i) = "alb_num_element";

i = i + 1;
kw_label(i) = "penalty_alpha";

i = i + 1;
kw_label(i) = "element_size";

i = i + 1;
kw_label(i) = "dg_solver";


i = i + 1;
kw_label(i) = "ecut_wavefunction";

i = i + 1;
kw_label(i) = "density_grid_factor";

i = i + 1;
kw_label(i) = "lgl_grid_factor";

i = i + 1;
kw_label(i) = "distance_periodize";

i = i + 1;
kw_label(i) = "buffer_size";


end