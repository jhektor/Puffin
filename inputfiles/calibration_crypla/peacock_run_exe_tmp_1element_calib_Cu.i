# [Mesh]
# type = FileMesh
# file = box-4x4x2um.inp
# displacements = 'disp_x disp_y disp_z'
# uniform_refine = 0
# []
# [Variables]
# [./disp_x]
# #scaling = 1e-2
# [../]
# [./disp_y]
# #scaling = 1e-2
# [../]
# [./disp_z]
# #scaling = 1e-2
# [../]
# []
# Tensor Mechanics action, sets up kernels, strain material and variables
[Mesh]
  type = GeneratedMesh
  displacements = 'disp_x disp_y disp_z'
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
  zmax = 50 # 8
[]

[GlobalParams]
  # CahnHilliard needs the third derivatives
  # enable_jit = true
  displacements = 'disp_x disp_y disp_z'
[]

[BCs]
  [symmy]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [symmx]
    type = PresetBC
    variable = disp_x
    boundary = 'left' # left right'
    value = 0
  []
  [symmz]
    type = PresetBC
    variable = disp_z
    boundary = 'back' # front back'
    value = 0
  []
  [loading]
    # function = '-0.005*t' #should give 1e-4 strain rate
    type = FunctionPresetBC
    variable = disp_z
    boundary = 'front'
    function = -0.05*t # should give 1e-3 strain rate as in Y. Kariya et. al. 2012
  []
[]

[UserObjects]
  # Crystal plasticity for central Sn grain
  [slip_rate_gss]
    # variable_size =10# 32
    # variable_size = 32
    # slip_sys_file_name = slip_systems_bct_v2_sorted.txt
    # flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1 13 32 1 1' #start_ss end_ss gamma0 1/m
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 12
    slip_sys_file_name = slip_sys_FCC.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 12 0.001 0.012' # start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss
  []
  [slip_resistance_gss]
    # variable_size = 10 #32
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  []
  [state_var_gss]
    # groups = '0 3 6 9 12'
    # group_values = '0.144 0.07 0.01' # 23 MPa in eV/nm^3 initial values of slip resistance
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'
    group_values = '0.1'
    uo_state_var_evol_rate_comp_name = 'state_var_evol_rate_comp_gss'
    scale_factor = '1.0'
  []
  [state_var_evol_rate_comp_gss]
    # groups = '0 4 8 12 16 20 32'
    # hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    # hprops = '1.4 0.624 0.250 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    # tau0_group_values = '0.12 0.12 0.12 0.12 0.12 0.12' # shold be set to zero as default -> Kalidindi's
    # tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0' # shold be set to zero as default -> Kalidindi's formulation
    type = CrystalPlasticityStateVarRateComponentVoce
    variable_size = 12
    groups = '0 12'
    h0_group_values = '1.125'
    tau0_group_values = '0.0'
    tauSat_group_values = '0.925'
    hardeningExponent_group_values = '2.25'
    coplanarHardening_group_values = '1.0' # q_aa = 1
    selfHardening_group_values = '1.4'
    crystal_lattice_type = FCC
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  []
[]

[Materials]
  # central Sn
  # [./strain]
  # type = ComputeFiniteStrain
  # #displacements = 'disp_x disp_y disp_z'
  # [../]
  [crysp]
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-3
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    maximum_substep_iteration = 20
    tan_mod_type = exact
  []
  [elasticity_tensor]
    # C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table I
    # C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 149.8022 137.31 137.31' #eV/nm^3 #From Darbandi 2013 table I
    # angles are ZX'Z'' rotations, proper euler, according to http://mooseframework.org/docs/doxygen/modules/RotationTensor_8C_source.html
    type = ComputeElasticityTensorCPBaseName # Allows for changes due to crystal re-orientation
    C_ijkl = '116100.0 58050.0 58050.0 116100.0 58050.0 116100.0 29030.0 29030.0 29030.0' # eV/nm^3 #From Kalidindi 1992
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  []
[]

[Debug]
  # show_material_props = true
  show_var_residual_norms = true
[]

[Executioner]
  # end_time = 1200 #Should give 12% strain for Darbandi 2013?
  # dt = 12
  type = Transient
  end_time = 15 # Should give 5% strain for 5comp  dir Y. Kariya 2012
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
  nl_max_its = 20
  [TimeStepper]
    # Turn on time stepping
    # postprocessor_dtlim = 5
    type = IterationAdaptiveDT
    dt = 0.1
    cutback_factor = 0.5
    growth_factor = 1.5
    optimal_iterations = 10
    linear_iteration_ratio = 10
  []
[]

[Outputs]
  file_base = calibrationSn
  exodus = true
  csv = true
  console = true
  print_linear_residuals = false
  print_perf_log = false
[]

[AuxVariables]
  [slip0]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip1]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip2]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip3]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip4]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip5]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip6]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip7]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip8]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip9]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip10]
    family = MONOMIAL
    order = CONSTANT
  []
  [slip11]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [slip0]
    type = MaterialStdVectorAux
    variable = slip0
    index = 0
    property = slip_rate_gss
  []
  [slip1]
    type = MaterialStdVectorAux
    variable = slip1
    index = 1
    property = slip_rate_gss
  []
  [slip2]
    type = MaterialStdVectorAux
    variable = slip2
    index = 2
    property = slip_rate_gss
  []
  [slip3]
    type = MaterialStdVectorAux
    variable = slip3
    index = 3
    property = slip_rate_gss
  []
  [slip4]
    type = MaterialStdVectorAux
    variable = slip4
    index = 4
    property = slip_rate_gss
  []
  [slip5]
    type = MaterialStdVectorAux
    variable = slip5
    index = 5
    property = slip_rate_gss
  []
  [slip6]
    type = MaterialStdVectorAux
    variable = slip6
    index = 6
    property = slip_rate_gss
  []
  [slip7]
    type = MaterialStdVectorAux
    variable = slip7
    index = 7
    property = slip_rate_gss
  []
  [slip8]
    type = MaterialStdVectorAux
    variable = slip8
    index = 8
    property = slip_rate_gss
  []
  [slip9]
    type = MaterialStdVectorAux
    variable = slip9
    index = 9
    property = slip_rate_gss
  []
  [slip10]
    type = MaterialStdVectorAux
    variable = slip10
    index = 10
    property = slip_rate_gss
  []
  [slip11]
    type = MaterialStdVectorAux
    variable = slip11
    index = 11
    property = slip_rate_gss
  []
[]

[Postprocessors]
  [slip0]
    type = PointValue
    variable = slip0
    point = '25 25 25'
  []
  [slip1]
    type = PointValue
    variable = slip1
    point = '25 25 25'
  []
  [slip2]
    type = PointValue
    variable = slip2
    point = '25 25 25'
  []
  [slip3]
    type = PointValue
    variable = slip3
    point = '25 25 25'
  []
  [slip4]
    type = PointValue
    variable = slip4
    point = '25 25 25'
  []
  [slip5]
    type = PointValue
    variable = slip5
    point = '25 25 25'
  []
  [slip6]
    type = PointValue
    variable = slip6
    point = '25 25 25'
  []
  [slip7]
    type = PointValue
    variable = slip7
    point = '25 25 25'
  []
  [slip8]
    type = PointValue
    variable = slip8
    point = '25 25 25'
  []
  [slip9]
    type = PointValue
    variable = slip9
    point = '25 25 25'
  []
  [slip10]
    type = PointValue
    variable = slip10
    point = '25 25 25'
  []
  [slip11]
    type = PointValue
    variable = slip11
    point = '25 25 25'
  []
  [stress_center_zz]
    type = PointValue
    variable = stress_zz
    point = '25 25 25'
  []
  [strain_center_zz]
    type = PointValue
    variable = strain_zz
    point = '25 25 25'
  []
  [stress_center_xy]
    type = PointValue
    variable = stress_xy
    point = '25 25 25'
  []
  [strain_center_xy]
    type = PointValue
    variable = strain_xy
    point = '25 25 25'
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [block1]
        strain = FINITE
        add_variables = true
        generate_output = 'stress_zz stress_xy strain_zz strain_xy vonmises_stress'
      []
    []
  []
[]
