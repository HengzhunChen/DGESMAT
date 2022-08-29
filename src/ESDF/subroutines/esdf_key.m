function esdf_key()
% ESDF_KEY Module to hold keyword list, used to initialize global variable
%    kw_label.
%
%    This must be updated as new keywords are brought into existence.
%
%    See also esdf_init, ESDFReadInput, ESDFInputParam.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
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
kw_label(i) = "restart_mode";

i = i + 1;
kw_label(i) = "restart_density";

i = i + 1;
kw_label(i) = "restart_wfn";


i = i + 1;
kw_label(i) = "output_dir";

i = i + 1;
kw_label(i) = "output_density";

i = i + 1;
kw_label(i) = "output_potential";

i = i + 1;
kw_label(i) = "output_wfn";

i = i + 1;
kw_label(i) = "output_alb_elem_lgl";

i = i + 1;
kw_label(i) = "output_alb_elem_uniform";

i = i + 1;
kw_label(i) = "output_wfn_extelem";

i = i + 1;
kw_label(i) = "output_pot_extelem";

i = i + 1;
kw_label(i) = "output_hmatrix";

i = i + 1;
kw_label(i) = "output_struct_info";


i = i + 1;
kw_label(i) = "potential_barrier";

i = i + 1;
kw_label(i) = "periodize_potential";

i = i + 1;
kw_label(i) = "calculate_aposteriori_each_scf";

i = i + 1;
kw_label(i) = "calculate_force_each_scf";


i = i + 1;
kw_label(i) = "restart_position";

i = i + 1;
kw_label(i) = "restart_velocity";

i = i + 1;
kw_label(i) = "output_position";

i = i + 1;
kw_label(i) = "output_velocity";

i = i + 1;
kw_label(i) = "output_xyz";



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
kw_label(i) = "ion_temperature";

i = i + 1;
kw_label(i) = "extra_electron";


i = i + 1;
kw_label(i) = "statfile";



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
kw_label(i) = "element_position_start";

i = i + 1;
kw_label(i) = "element_grid_size";

i = i + 1;
kw_label(i) = "dg_solver";


i = i + 1;
kw_label(i) = "potential_barrier_w";

i = i + 1;
kw_label(i) = "potential_barrier_s";

i = i + 1;
kw_label(i) = "potential_barrier_r";

i = i + 1;
kw_label(i) = "ecut_wavefunction";

i = i + 1;
kw_label(i) = "density_grid_factor";

i = i + 1;
kw_label(i) = "lgl_grid_factor";

i = i + 1;
kw_label(i) = "gauss_interp_factor";

i = i + 1;
kw_label(i) = "gauss_sigma";

i = i + 1;
kw_label(i) = "distance_periodize";



% Inputs related to Chebyshev ploynomial Filtered SCF iterations for DG
i = i + 1;
kw_label(i) = "first_scfdg_chebyfilterorder";

i = i + 1;
kw_label(i) = "first_scfdg_chebycyclenum";

i = i + 1;
kw_label(i) = "second_scfdg_chebyouteriter";

i = i + 1;
kw_label(i) = "second_scfdg_chebyfilterorder";

i = i + 1;
kw_label(i) = "second_scfdg_chebycyclenum";

i = i + 1;
kw_label(i) = "general_scfdg_chebyfilterorder";

i = i + 1;
kw_label(i) = "general_scfdg_chebycyclenum";


% Inputs related to Chebyshev polynomial filtered 
% complementary subspace iteration strategy in DGDFT
i = i + 1;
kw_label(i) = "scfdg_use_chefsi_complementary_subspace";

i = i + 1;
kw_label(i) = "scfdg_chefsi_complementary_subspace_syrk";

i = i + 1;
kw_label(i) = "scfdg_chefsi_complementary_subspace_syr2k";

i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_nstates";

i = i + 1;
kw_label(i) = "scfdg_cs_ioniter_regular_cheby_freq";

i = i + 1;
kw_label(i) = "scfdg_cs_bigger_grid_dim_fac";


% Inner LOBPCG related options
i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_inner_lobpcgtol";

i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_inner_lobpcgiter";


% Inner CheFSI related options
i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_use_inner_cheby";

i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_inner_chebyfilterorder";

i = i + 1;
kw_label(i) = "scfdg_complementary_subspace_inner_chebycyclenum";




% ************************************************************************
%                    Parameters for Ionic Motion
% ************************************************************************

i = i + 1;
kw_label(i) = "ion_energy_diff";

i = i + 1;
kw_label(i) = "md_scf_energy_criteria_engage_ioniter";

i = i + 1;
kw_label(i) = "md_scf_outer_maxiter";

i = i + 1;
kw_label(i) = "md_scf_etot_diff";

i = i + 1;
kw_label(i) = "md_scf_eband_diff";


i = i + 1;
kw_label(i) = "energy_gap";

i = i + 1;
kw_label(i) = "spectral_radius";

i = i + 1;
kw_label(i) = "matrix_ordering";

i = i + 1;
kw_label(i) = "inertia_count";

i = i + 1;
kw_label(i) = "inertia_count_steps";

i = i + 1;
kw_label(i) = "unused_states";


i = i + 1;
kw_label(i) = "ion_max_iter";

i = i + 1;
kw_label(i) = "ion_move";

i = i + 1;
kw_label(i) = "geo_opt_max_force";

i = i + 1;
kw_label(i) = "geo_opt_nlcg_sigma";


i = i + 1;
kw_label(i) = "fire_nmin";

i = i + 1;
kw_label(i) = "fire_time_step";

i = i + 1;
kw_label(i) = "fire_atomic_mass";

i = i + 1;
kw_label(i) = "md_max_step";

i = i + 1;
kw_label(i) = "md_time_step";

i = i + 1;
kw_label(i) = "md_extrapolation_type";

i = i + 1;
kw_label(i) = "md_extrapolation_variable";

i = i + 1;
kw_label(i) = "md_extrapolation_wavefunction";

i = i + 1;
kw_label(i) = "thermostat_mass";

i = i + 1;
kw_label(i) = "langevin_damping";

i = i + 1;
kw_label(i) = "kappa_xlbomd";


% ************************************************************************
%                        Parameters for Hybird
% ************************************************************************

i = i + 1;
kw_label(i) = "scf_phi_maxiter";

i = i + 1;
kw_label(i) = "md_scf_phi_maxiter";

i = i + 1;
kw_label(i) = "scf_phi_tolerance";

i = i + 1;
kw_label(i) = "hybrid_ace";

i = i + 1;
kw_label(i) = "hybrid_df";

i = i + 1;
kw_label(i) = "hybrid_df_type";

i = i + 1;
kw_label(i) = "hybrid_df_kmeans_wf_type";

i = i + 1;
kw_label(i) = "hybrid_df_kmeans_wf_alpha";

i = i + 1;
kw_label(i) = "hybrid_df_kmeans_tolerance";

i = i + 1;
kw_label(i) = "hybrid_df_kmeans_maxiter";

i = i + 1;
kw_label(i) = "hybrid_df_num_mu";

i = i + 1;
kw_label(i) = "hybrid_df_num_gaussianrandom";

i = i + 1;
kw_label(i) = "hybrid_df_tolerance";

i = i + 1;
kw_label(i) = "hybrid_active_init";

i = i + 1;
kw_label(i) = "hybrid_mixing_type";

i = i + 1;
kw_label(i) = "hybrid_ace_twice_pcdiis";

i = i + 1;
kw_label(i) = "exx_divergence_type";



end