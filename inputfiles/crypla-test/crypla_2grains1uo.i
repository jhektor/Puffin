[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 10
  ny = 10
  nz = 1
  xmin = -0.5
  xmax = 0.5
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

  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./vm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.01*t
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
    function = '1-0.5*(1+tanh(x/0.01))'
    use_displaced_mesh = true
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
    scalar_type = VonMisesStress
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = syy
    rank_two_tensor = global_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./stress_yyf]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = global_stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = global_total_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = global_fp
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
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = disp_x
    boundary = bottom
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = bottom
    value = 0
  [../]
  [./tdisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = tdisp
  [../]
[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
    flowprops = '1 32 0.001 0.05' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss
    base_name = 'global'
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    groups = '0 32'
    group_values = '23' # MPa initial values of slip resistance
    #group_values = '23000' # MPa to make everything elastic
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 32
    hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]

[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCPBaseName
    block = 0
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    base_name = 'global'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
    base_name = 'global'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    block = 0
    C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    base_name = 'global'
  [../]

[]

[Postprocessors]
  [./e_yy1]
    type = PointValue
    variable = e_yy
    point = '0 0.5 0.05'
    #use_displaced_mesh = true
  [../]
  [./stress_yy1]
    type = PointValue
    variable = stress_yy
    point = '0 0.5 0.05'
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

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  nl_abs_tol = 1e-10
  nl_rel_step_tol = 1e-10
  dtmax = 10.0
  nl_rel_tol = 1e-10
  ss_check_tol = 1e-10
  #end_time = 1
  #dtmin = 0.05
  num_steps = 10000
  nl_abs_step_tol = 1e-10
  l_tol = 1e-4
[]

[Outputs]
  file_base = 2grains-elastic-tension
  exodus = true
[]
