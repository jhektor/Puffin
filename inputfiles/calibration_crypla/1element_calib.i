[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 1
  ny = 1
  nz = 1
  xmin = 0
  xmax = 50
  ymin = 0
  ymax = 50
  zmin = 0
  zmax = 50 #8
  displacements = 'disp_x disp_y disp_z'
[]
#[Mesh]
#  type = FileMesh
#  file = box-4x4x2um.inp
#  displacements = 'disp_x disp_y disp_z'
#  uniform_refine = 0
#[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  #enable_jit = true
  displacements = 'disp_x disp_y disp_z'
[]
#[Variables]
#  [./disp_x]
#    #scaling = 1e-2
#  [../]
#  [./disp_y]
#    #scaling = 1e-2
#  [../]
#  [./disp_z]
#    #scaling = 1e-2
#  [../]
#[]

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
    boundary = 'left'# left right'
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = 'back'# front back'
    value = 0
  [../]
  [./loading]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = '0.007*t'
  [../]

[]
[UserObjects]
  #Crystal plasticity for central Sn grain
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
    flowprops = '1 32 0.001 0.05' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss
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
    group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 32
    #hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    hprops = '1.4 0.624 0.250 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
  []

[Materials]
  #central Sn
  [./crysp]
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    maximum_substep_iteration = 10
    tan_mod_type = exact
  [../]
  #[./strain]
  #  type = ComputeFiniteStrain
  #  #displacements = 'disp_x disp_y disp_z'
  #[../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 149.8022 137.31 137.31' #eV/nm^3 #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 60
    euler_angle_2 = 90
    euler_angle_3 = 0
  [../]

[]

#Tensor Mechanics action, sets up kernels, strain material and variables
[Modules/TensorMechanics/Master]
  [./block1]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_yy stress_xz strain_yy strain_xz vonmises_stress'
  [../]
[]


[Postprocessors]
  [./stress_center_yy]
    type = PointValue
    variable = stress_yy
    point = '25 25 25'
  [../]
  [./strain_center_yy]
    type = PointValue
    variable = strain_yy
    point = '25 25 25'
  [../]
  [./stress_center_xz]
    type = PointValue
    variable = stress_xz
    point = '25 25 25'
  [../]
  [./strain_center_xz]
    type = PointValue
    variable = strain_xz
    point = '25 25 25'
  [../]
[]
[Debug]
  show_var_residual_norms = true
  #show_material_props = true
[]

[Executioner]
  type = Transient
  end_time = 16
  dt = 0.25
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
[]

[Outputs]
  file_base = calibrationSn
  exodus = true
  csv = true
  print_perf_log = true
[]
