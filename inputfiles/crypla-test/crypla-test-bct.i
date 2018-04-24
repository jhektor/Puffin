[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 3
  ny = 3
  nz = 3
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
  [./eyy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./sxx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./syy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./szz]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.001*t
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./eyy]
    type = RankTwoAux
    rank_two_tensor = lage
    variable = eyy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./sxx]
    type = RankTwoAux
    rank_two_tensor = pk2
    variable = sxx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]
  [./syy]
    type = RankTwoAux
    rank_two_tensor = pk2
    variable = syy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./szz]
    type = RankTwoAux
    rank_two_tensor = pk2
    variable = szz
    index_i = 2
    index_j = 2
    execute_on = 'initial timestep_end'
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
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = back
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
  #[./slip_rate_gss]
  #  type = CrystalPlasticitySlipRateGSS
  #  variable_size = 32
  #  slip_sys_file_name = slip_systems_bct.txt
  #  num_slip_sys_flowrate_props = 2
  #  #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
  #  flowprops = '1 32 0.001 0.05' #start_ss end_ss gamma0 1/m
  #  uo_state_var_name = state_var_gss
  #[../]
  #[./slip_resistance_gss]
  #  type = CrystalPlasticitySlipResistanceGSS
  #  variable_size = 32
  #  uo_state_var_name = state_var_gss
  #[../]
  #[./state_var_gss]
  #  type = CrystalPlasticityStateVariable
  #  variable_size = 32
  #  groups = '0 32'
  #  #group_values = '23' # MPa initial values of slip resistance
  #  group_values = '23000' # MPa initial values of slip resistance
  #  uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
  #  scale_factor = 1.0
  #[../]
  #[./state_var_evol_rate_comp_gss]
  #  type = CrystalPlasticityStateVarRateComponentGSS
  #  variable_size = 32
  #  hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
  #  uo_slip_rate_name = slip_rate_gss
  #  uo_state_var_name = state_var_gss
  #[../]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSS
    variable_size = 12
    slip_sys_file_name = input_slip_sys.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    uo_state_var_name = state_var_gss
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 4 8 12'
    group_values = '60000.8 60000.8 60000.8'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1.0 541.5 109.8 2.5'
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCP
    block = 0
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP #Allows for changes due to crystal re-orientation
    block = 0
    C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2 = 90
    euler_angle_3 = 90
  [../]
[]

[Postprocessors]
  [./eyy]
    type = ElementAverageValue
    variable = eyy
    point = '0.5 0.5 0.5'
  [../]
  [./sxx]
    type = ElementAverageValue
    variable = sxx
    point = '0.5 0.5 0.5'
  [../]
  [./syy]
    type = ElementAverageValue
    variable = syy
    point = '0.5 0.5 0.5'
  [../]
  [./szz]
    type = ElementAverageValue
    variable = szz
    point = '0.5 0.5 0.5'
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
  num_steps = 200
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  file_base = tension-45-90-90-elastic
  exodus = true
  csv = true
[]
