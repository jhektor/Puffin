[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 10
  ny = 10
  nz = 1
  xmin = -800
  xmax = 800
  ymin = -800
  ymax = 800
  zmin = 0
  zmax = 16
  displacements = 'disp_x disp_y disp_z'
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  displacements = 'disp_x disp_y disp_z'
[]
[Variables]
  [./disp_x]
    scaling = 1e-2
  [../]
  [./disp_y]
    scaling = 1e-2
  [../]
  [./disp_z]
    scaling = 1e-2
  [../]
[]
[BCs]
  [./Periodic]
    #generated mesh
    [./xy]
      auto_direction = 'x z'
      variable = 'disp_x disp_y disp_z'
    [../]
  [../]
  [./symmy]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = disp_x
    boundary = 'bottom'# left right'
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = 'bottom'# front back'
    value = 0
  [../]
  [./load]
    type = FunctionPresetBC
    boundary = 'top'
    variable = 'disp_y'
    function = '-t'
  [../]
[]
[AuxVariables]
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./shyd]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
[AuxKernels]
  [./eta]
    type = FunctionAux
    variable = eta
    function = 'if(x<0,1,0)'
  [../]
  [./eta2]
    type = FunctionAux
    variable = eta2
    function = 'if(x>=0,1,0)'
  [../]
  [./shyd]
    type = RankTwoScalarAux
    rank_two_tensor = global_stress
    variable = shyd
    scalar_type = Hydrostatic
  [../]
[]
[UserObjects]
  #Crystal plasticity for central Sn grain
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct_miller.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    #flowprops = '1 32 0.001 2.6e-8' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss
    base_name = 'eta'
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    #groups = '0 32'
    #group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentVoceP
    variable_size = 32
    groups = '0 2 4 6 10 20 24 32'
    h0_group_values = '0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484'
    tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    tauSat_group_values = '0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115'
    hardeningExponent_group_values = '2.0 2.0 2.0 2.0 2.0 2.0 2.0'
    coplanarHardening_group_values = '1.0 1.0 1.0 1.0 1.0 1.0 1.0' #q_aa = 1
    selfHardening_group_values = '1.4 1.4 1.4 1.4 1.4 1.4 1.4'

    crystal_lattice_type = BCT # default is BCT
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
  [./slip_rate_gss2]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct_miller.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    #flowprops = '1 32 0.001 2.6e-8' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss2
    base_name = 'eta2'
  [../]
  [./slip_resistance_gss2]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss2
  [../]
  [./state_var_gss2]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    #groups = '0 32'
    #group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss2
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss2]
    type = CrystalPlasticityStateVarRateComponentVoceP
    variable_size = 32
    groups = '0 2 4 6 10 20 24 32'
    h0_group_values = '0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484'
    tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    tauSat_group_values = '0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115'
    hardeningExponent_group_values = '2.0 2.0 2.0 2.0 2.0 2.0 2.0'
    coplanarHardening_group_values = '1.0 1.0 1.0 1.0 1.0 1.0 1.0' #q_aa = 1
    selfHardening_group_values = '1.4 1.4 1.4 1.4 1.4 1.4 1.4'

    crystal_lattice_type = BCT # default is BCT
    uo_slip_rate_name = slip_rate_gss2
    uo_state_var_name = state_var_gss2
  [../]
[]
[Materials]
  [./h]
    type = SwitchingFunctionMaterial
    eta = eta
    output_properties = h
    outputs = exodus
    function_name = h
  [../]
  [./h2]
    type = SwitchingFunctionMaterial
    eta = eta2
    output_properties = h2
    outputs = exodus
    function_name = h2
  [../]
  [./crysp]
    # type = FiniteStrainUObasedCPBaseNamehscale
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    base_name = 'eta'
    maximum_substep_iteration = 10
    tan_mod_type = exact
    h_scale = 'h'
    h_tol = 0.01
  [../]
  [./strain]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = 'eta'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3' #MPa #From Darbandi 2014 table I
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 90
    euler_angle_2 = 15
    euler_angle_3 = 0
    base_name = 'eta'
  [../]

  [./crysp2]
    # type = FiniteStrainUObasedCPBaseNamehscale
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate_gss2'
    uo_slip_resistances = 'slip_resistance_gss2'
    uo_state_vars = 'state_var_gss2'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss2'
    base_name = 'eta2'
    maximum_substep_iteration = 10
    tan_mod_type = exact
    h_scale = 'h2'
    h_tol = 0.01
  [../]
  [./strain2]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = 'eta2'
  [../]
  [./elasticity_tensor2]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3' #MPa #From Darbandi 2014 table I
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 90
    euler_angle_2 = 90
    euler_angle_3 = 90
    base_name = 'eta2'
  [../]


  #
  # [./elasticity_tensor2]
  #   type = ComputeIsotropicElasticityTensor
  #   youngs_modulus = 400. #150 GPa in eV/nm^3
  #   poissons_ratio = 0.35
  #   base_name = eta2
  # [../]
  # [./strain2]
  #   type = ComputeFiniteStrain
  #   base_name = eta2
  # [../]
  # [./stress2]
  #   type = ComputeFiniteStrainElasticStress
  #   base_name = eta2
  # [../]

  [./global_stress] #homogeniserar bara Cauchy stress
    type = MultiPhaseStressMaterial
    phase_base = 'eta eta2'
    h          = 'h h2'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = global
  [../]
[]
[Kernels]
  #Stress divergence
  [./TensorMechanics]
    #displacements = 'disp_x disp_y disp_z'
    base_name = global
    use_displaced_mesh = true
    strain = FINITE #this also sets incremental strain =true
    generate_output = 'stress_xx stress_yy stress_zz'
  [../]
[]

[AuxVariables]
  [./slips]
    order = FIRST
    family = MONOMIAL
  [../]
  [./accum_slip]
    order = FIRST
    family = MONOMIAL
  [../]
[]
[AuxKernels]
  [./slips]
    type =  CrystalPlasticityTotalSlip
    variable = slips
    h_names = 'h2'
    slip_rates = 'slip_rate_gss2'
    execute_on = 'initial TIMESTEP_END'
  [../]
  [./accumslip]
    type = AccumulateAux
    accumulate_from_variable = slips
    variable = accum_slip
    execute_on = 'initial TIMESTEP_END'
  [../]
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
[Preconditioning]
  [./smp]
    type = SMP
    full = true
    solve_type = PJFNK

    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu       ilu             nonzero             nonzero            superlu_dist'
    #

    #petsc_options_iname = '-pc_asm_overlap -ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '1  121              1e-10     1e-8     1e-5          asm       ilu             nonzero             nonzero'
    #petsc_options_iname = '-pc_asm_overlap -ksp_gmres_restart  -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     1               121                asm       ilu             nonzero             nonzero'

    #Larrys suggestion
    petsc_options_iname = '-pc_asm_overlap -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type -sub_ksp_type'
    petsc_options_value = '2                  asm       lu             nonzero             nonzero                    preonly'

    #petsc_options_iname = '-pc_type -pc_factor_shift_type -sub_pc_factor_shift_type'
    #petsc_options_value = 'bjacobi       nonzero  nonzero'

    #petsc_options_iname = '-pc_type -pc_hypre_type   -pc_factor_shift_type  -ksp_gmres_restart'
    #petsc_options_value = 'hypre       boomeramg            nonzero    150'
  [../]
[]

[Executioner]
   type = Transient
   solve_type = 'PJFNK'
   #solve_type = 'NEWTON'

   line_search = basic
   #line_search = none
   #line_search = default
   #line_search = bt

   petsc_options = '-snes_converged_reason -ksp_converged_reason
   -snes_ksp_ew'

   #l_max_its = 120
   l_max_its = 10
   nl_max_its = 25
   #nl_max_its = 15

   l_tol = 1.0e-3
   #l_abs_step_tol = 1e-8
   nl_rel_tol = 1.0e-6 #1.0e-9 # 1.0e-8 # 1.0e-9 #1.0e-10
   nl_abs_tol = 1.0e-9 #1.0e-10 # 1.0e-9 # 1.0e-10#1.0e-11

   #num_steps = 2000
   end_time = 4
   # end_time = 100 #2x176h
   #n_startup_steps = 2
   #dtmin= 1e-5
   dtmax= 1e5

   #scheme = implicit-euler
   #scheme = bdf2
   # scheme = dirk
   [./TimeIntegrator]
     # type = ExplicitTVDRK2
       type = LStableDirk2
     # type = ImplicitEuler
   [../]

   [./TimeStepper]
       # Turn on time stepping
       type = IterationAdaptiveDT
       #dt = 1e-3
       dt = 5e-3
       #dt = 2e-2


       cutback_factor = 0.5
       growth_factor = 1.3
       #optimal_iterations = 10
       optimal_iterations = 15 # 20
       linear_iteration_ratio = 25
       #postprocessor_dtlim = 5
   [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = scaletest #1Cu1imc2sn-initial_IMC_001
    append_date = true
  [../]
  [./perf_log]
    type =  CSV
    file_base = perf_log_scaled
    append_date = true
  [../]
  [./sim_log]
    type = Console
    append_date = true
    output_file = true
    file_base = console_sim_log_scaled
  [../]
  perf_graph = true
[]
