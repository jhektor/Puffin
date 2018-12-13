[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 30
  ny = 30
  nz = 3
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = 0
  zmax = 0.1
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    block = 0
  [../]
  [./disp_y]
    block = 0
  [../]
  [./disp_z]
    block = 0
  [../]
[]

[AuxVariables]
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]

  [./e_yy1]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_yy1]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_yy2]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_yy2]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]

  [./vm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syy1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syy2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syyg]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxyg]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxzg]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eyyg]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ep_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.01*t
  [../]
  #[./teta]
  #  type = ParsedFunction
  #  value = 1-0.5*(1+tanh((x)/0.1))
  #[../]
  [./teta]
    type = ParsedFunction
    #value = 'r:=sqrt(x^2+y^2);if(r<=0.1*t,0,1)'
    value = 'r:=sqrt(x^2+y^2);0.5*(1+tanh((r-0.1*t)/0.1))'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    base_name = global
  [../]
[]

[AuxKernels]
  [./bnd1]
    type = FunctionAux
    variable = eta1
    function = teta
    use_displaced_mesh = false
  [../]
  [./bnd2]
    type = ParsedAux
    variable = eta2
    args = eta1
    function = 1-eta1
  [../]
  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = vm
    rank_two_tensor =global_stress
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]
  [./fp_xx1]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = eta1_fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yy1]
    type = RankTwoAux
    variable = stress_yy1
    rank_two_tensor = eta1_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yy2]
    type = RankTwoAux
    variable = stress_yy2
    rank_two_tensor = eta2_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yyg]
    type = RankTwoAux
    variable = syyg
    rank_two_tensor = global_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_xyg]
    type = RankTwoAux
    variable = sxyg
    rank_two_tensor = global_stress
    index_j = 0
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_xzg]
    type = RankTwoAux
    variable = sxzg
    rank_two_tensor = global_stress
    index_j = 0
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yy1f]
    type = RankTwoAux
    variable = syy1
    rank_two_tensor = eta1_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yy2f]
    type = RankTwoAux
    variable = syy2
    rank_two_tensor = eta2_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yy1]
    type = RankTwoAux
    variable = e_yy1
    rank_two_tensor = eta1_total_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yy2]
    type = RankTwoAux
    variable = e_yy2
    rank_two_tensor = eta2_total_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yyg]
    type = RankTwoAux
    variable = eyyg
    rank_two_tensor = global_total_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./symmy]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = disp_x
    boundary = 'left right'
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = 'front back'
    value = 0
  [../]
  #[./tdisp]
  #  type = FunctionPresetBC
  #  variable = disp_y
  #  boundary = top
  #  function = tdisp
  #[../]
[]

[UserObjects]
  [./slip_rate_gss1]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
    flowprops = '1 32 0.001 0.05' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss1
    base_name = 'eta1'
  [../]
  [./slip_resistance_gss1]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss1
  [../]
  [./state_var_gss1]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    groups = '0 32'
    group_values = '23' # MPa initial values of slip resistance
    #group_values = '23000' # MPa to make everything elastic
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss1
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss1]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 32
    hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    uo_slip_rate_name = slip_rate_gss1
    uo_state_var_name = state_var_gss1
  [../]
[]

[Materials]
  [./time]
    type = TimeStepMaterial
    prop_time = time
    prop_dt = dt
    use_displaced_mesh = true
  [../]
  [./crysp]
    type = FiniteStrainUObasedCPBaseName
    block = 0
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss1'
    uo_slip_resistances = 'slip_resistance_gss1'
    uo_state_vars = 'state_var_gss1'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss1'
    base_name = 'eta1'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
    base_name = 'eta1'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    block = 0
    C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    base_name = 'eta1'
  [../]
  [./crysp2]
    type = ComputeFiniteStrainElasticStress
    base_name = 'eta2'
  [../]
  [./strain2]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
    base_name = 'eta2'
    eigenstrain_names = eT_eps
  [../]
  [./elasticity_tensor2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 112300. #112.3 GPa
    poissons_ratio = 0.31
    base_name = eta2
  [../]
  [./eigenstrain_eps]
    type = ComputeVariableEigenstrain
    #args = 'y'
    args = 'eta1 eta2'
    base_name = eta2
    eigen_base = '1 1 1 0 0 0'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eps
    prefactor = pre_eps
  [../]
  [./pre_eps]
    type = DerivativeParsedMaterial
    f_name = pre_eps
    material_property_names = 'h2 time'
    function = '0.02*h2*if(time>1,1,time)'
    outputs = exodus
    output_properties = pre_eps
  [../]
  #[./pre_eps]
  #  type = GenericFunctionMaterial
  #  prop_names = pre_eps
  #  prop_values = tdisp
  #  #prop_values = 0.02
  #  outputs = exodus
  #  output_properties = pre_eps
  #[../]
  [./h1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta1 eta2'
      phase_etas = eta1
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h1
  [../]
  [./h2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta1 eta2'
      phase_etas = eta2
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h2
  [../]
  # Generate the global stress from the phase stresses
  [./global_stress] #homogeniserar bara Cauchy stress
    type = MultiPhaseStressMaterial
    phase_base = 'eta1 eta2'
    h          = 'h1 h2'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
    base_name = global
  [../]
[]

[Postprocessors]
  [./e_yy1]
    type = PointValue
    variable = e_yy1
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_yy1]
    type = PointValue
    variable = stress_yy1
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./e_yy2]
    type = PointValue
    variable = e_yy2
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_yy2]
    type = PointValue
    variable = stress_yy2
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./e_yyg]
    type = PointValue
    variable = eyyg
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_yyg]
    type = PointValue
    variable = syyg
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_xyg]
    type = PointValue
    variable = stress_yy1
    point = '-0.25 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_xzg]
    type = PointValue
    variable = sxzg
    point = '0.0 0.5 0.05'
    #use_displaced_mesh = true
  [../]

[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.05

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  #petsc_options_iname = -pc_hypre_type
  #petsc_options_value = boomerang

  petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
  petsc_options_value = '     121              5e-9     1e-8     1e-5          lu       ilu             nonzero             nonzero            superlu_dist'

  nl_abs_tol = 5e-9
  nl_rel_step_tol = 1e-10
  dtmax = 10.0
  nl_rel_tol = 1e-8
  ss_check_tol = 1e-10
  #end_time = 1
  #dtmin = 0.05
  num_steps = 10000
  #nl_abs_step_tol = 5e-9
  l_tol = 1e-4
  nl_max_its = 20
[]

[Outputs]
  file_base = 2phase-eigenstrain-sharpcircle
  exodus = true
[]
