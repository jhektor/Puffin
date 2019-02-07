[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 1
  ny = 1
  nz = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 1 #8
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
    variable = disp_x
    boundary = 'right'
    # function = '-0.05*t' #should give 1e-3 strain rate as in Y. Kariya et. al. 2012
    function = '0.1*t-0.00022'
  [../]

[]

[UserObjects]
  #Crystal plasticity for central Sn grain
  [./slip_rate]
    type = CrystalPlasticitySlipRateHWR
    variable_size = 32
    gamma0 = 0.001
    m = 6
    slip_sys_file_name = slip_systems_bct_HWR.txt
    uo_slip_res_name = slip_resistance
  [../]
  [./slip_resistance]
    type = CrystalPlasticitySlipResistanceHWR
    #variable_size = 10 #32
    variable_size = 32
    groups = '0 32'
    G0 = '10'
    q = 1.4
    Q = 10
    crystal_lattice_type = 2
    uo_state_var_name = state_var
  [../]
  [./state_var]
    type = CrystalPlasticityStateVariableHWR
    variable_size = 32
    groups = '0 32'
    group_values = '0'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate]
    type = CrystalPlasticityStateVarRateComponentHWR
    variable_size = 32
    groups = '0 32'
    B = '8'
    uo_slip_resistance_name = slip_resistance
    uo_slip_rate_name = slip_rate
    uo_state_var_name = state_var
  [../]
[]
[Materials]
  #central Sn
  [./crysp]
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-3
    uo_slip_rates = 'slip_rate'
    uo_slip_resistances = 'slip_resistance'
    uo_state_vars = 'state_var'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate'
    maximum_substep_iteration = 20
    tan_mod_type = none
    output_properties = pk2_xx
    outputs = exodus
  [../]
  #[./strain]
  #  type = ComputeFiniteStrain
  #  #displacements = 'disp_x disp_y disp_z'
  #[../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    # C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table I
    C_ijkl = '19e3 0.36'
    fill_method = symmetric_isotropic_E_nu
    euler_angle_1 = 90
    euler_angle_2 = 60
    euler_angle_3 = 30
    # angles are ZX'Z'' rotations, proper euler, according to http://mooseframework.org/docs/doxygen/modules/RotationTensor_8C_source.html
  [../]
[]

#Tensor Mechanics action, sets up kernels, strain material and variables
[Modules/TensorMechanics/Master]
  [./block1]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx strain_xy vonmises_stress'
  [../]
[]

[Debug]
  show_var_residual_norms = true
  #show_material_props = true
[]

[Executioner]
  type = Transient
  #end_time = 1200 #Should give 12% strain for Darbandi 2013?
  num_steps = 100
  dt = 5e-3
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
  nl_max_its = 20



  # [./TimeStepper]
  #     # Turn on time stepping
  #     type = IterationAdaptiveDT
  #     dt = 0.05
  #     cutback_factor = 1.0
  #     growth_factor = 1.0
  #     optimal_iterations = 10
  #     linear_iteration_ratio = 10
  #     #postprocessor_dtlim = 5
  # [../]
[]

[Outputs]
  file_base = hwr
  exodus = true
  csv = true
  console = true
  print_linear_residuals = false
  perf_graph = true
[]

[AuxVariables]
  [./slip0]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip1]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip2]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip3]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip4]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip5]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip6]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip7]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip8]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip9]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip10]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip11]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip12]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip13]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip14]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip15]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip16]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip17]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip18]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip19]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip20]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip21]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip22]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip23]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip24]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip25]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip26]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip27]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip28]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip29]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip30]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./slip31]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./u_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./pk2_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./pk2_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./pk2_zz]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]
[AuxKernels]
  [./slip0]
    type = MaterialStdVectorAux
    variable = slip0
    index = 0
    property = slip_rate
  [../]
  [./slip1]
    type = MaterialStdVectorAux
    variable = slip1
    index = 1
    property = slip_rate
  [../]
  [./slip2]
    type = MaterialStdVectorAux
    variable = slip2
    index = 2
    property = slip_rate
  [../]
  [./slip3]
    type = MaterialStdVectorAux
    variable = slip3
    index = 3
    property = slip_rate
  [../]
  [./slip4]
    type = MaterialStdVectorAux
    variable = slip4
    index = 4
    property = slip_rate
  [../]
  [./slip5]
    type = MaterialStdVectorAux
    variable = slip5
    index = 5
    property = slip_rate
  [../]
  [./slip6]
    type = MaterialStdVectorAux
    variable = slip6
    index = 6
    property = slip_rate
  [../]
  [./slip7]
    type = MaterialStdVectorAux
    variable = slip7
    index = 7
    property = slip_rate
  [../]
  [./slip8]
    type = MaterialStdVectorAux
    variable = slip8
    index = 8
    property = slip_rate
  [../]
  [./slip9]
    type = MaterialStdVectorAux
    variable = slip9
    index = 9
    property = slip_rate
  [../]
  [./slip10]
    type = MaterialStdVectorAux
    variable = slip10
    index = 10
    property = slip_rate
  [../]
  [./slip11]
    type = MaterialStdVectorAux
    variable = slip11
    index = 11
    property = slip_rate
  [../]
  [./slip12]
    type = MaterialStdVectorAux
    variable = slip12
    index = 12
    property = slip_rate
  [../]
  [./slip13]
    type = MaterialStdVectorAux
    variable = slip13
    index = 13
    property = slip_rate
  [../]
  [./slip14]
    type = MaterialStdVectorAux
    variable = slip14
    index = 14
    property = slip_rate
  [../]
  [./slip15]
    type = MaterialStdVectorAux
    variable = slip15
    index = 15
    property = slip_rate
  [../]
  [./slip16]
    type = MaterialStdVectorAux
    variable = slip16
    index = 16
    property = slip_rate
  [../]
  [./slip17]
    type = MaterialStdVectorAux
    variable = slip17
    index = 17
    property = slip_rate
  [../]
  [./slip18]
    type = MaterialStdVectorAux
    variable = slip18
    index = 18
    property = slip_rate
  [../]
  [./slip19]
    type = MaterialStdVectorAux
    variable = slip19
    index = 19
    property = slip_rate
  [../]
  [./slip20]
    type = MaterialStdVectorAux
    variable = slip20
    index = 20
    property = slip_rate
  [../]
  [./slip21]
    type = MaterialStdVectorAux
    variable = slip21
    index = 21
    property = slip_rate
  [../]
  [./slip22]
    type = MaterialStdVectorAux
    variable = slip22
    index = 22
    property = slip_rate
  [../]
  [./slip23]
    type = MaterialStdVectorAux
    variable = slip23
    index = 23
    property = slip_rate
  [../]
  [./slip24]
    type = MaterialStdVectorAux
    variable = slip24
    index = 24
    property = slip_rate
  [../]
  [./slip25]
    type = MaterialStdVectorAux
    variable = slip25
    index = 25
    property = slip_rate
  [../]
  [./slip26]
    type = MaterialStdVectorAux
    variable = slip26
    index = 26
    property = slip_rate
  [../]
  [./slip27]
    type = MaterialStdVectorAux
    variable = slip27
    index = 27
    property = slip_rate
  [../]
  [./slip28]
    type = MaterialStdVectorAux
    variable = slip28
    index = 28
    property = slip_rate
  [../]
  [./slip29]
    type = MaterialStdVectorAux
    variable = slip29
    index = 29
    property = slip_rate
  [../]
  [./slip30]
    type = MaterialStdVectorAux
    variable = slip30
    index = 30
    property = slip_rate
  [../]
  [./slip31]
    type = MaterialStdVectorAux
    variable = slip31
    index = 31
    property = slip_rate
  [../]
  [./pk2_xx]
    type = MaterialRankTwoTensorAux
    variable = pk2_xx
    i = 0
    j = 0
    property = pk2
  [../]
  [./pk2_yy]
    type = MaterialRankTwoTensorAux
    variable = pk2_yy
    i = 1
    j = 1
    property = pk2
  [../]
  [./pk2_zz]
    type = MaterialRankTwoTensorAux
    variable = pk2_zz
    i = 2
    j = 2
    property = pk2
  [../]
[]
[Postprocessors]
  # [./slip0]
  #   type = ElementAverageValue
  #   variable = slip0
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip1]
  #   type = ElementAverageValue
  #   variable = slip1
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip2]
  #   type = ElementAverageValue
  #   variable = slip2
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip3]
  #   type = ElementAverageValue
  #   variable = slip3
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip4]
  #   type = ElementAverageValue
  #   variable = slip4
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip5]
  #   type = ElementAverageValue
  #   variable = slip5
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip6]
  #   type = ElementAverageValue
  #   variable = slip6
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip7]
  #   type = ElementAverageValue
  #   variable = slip7
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip8]
  #   type = ElementAverageValue
  #   variable = slip8
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip9]
  #   type = ElementAverageValue
  #   variable = slip9
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip10]
  #   type = ElementAverageValue
  #   variable = slip10
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip11]
  #   type = ElementAverageValue
  #   variable = slip11
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip12]
  #   type = ElementAverageValue
  #   variable = slip12
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip13]
  #   type = ElementAverageValue
  #   variable = slip13
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip14]
  #   type = ElementAverageValue
  #   variable = slip14
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip15]
  #   type = ElementAverageValue
  #   variable = slip15
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip16]
  #   type = ElementAverageValue
  #   variable = slip16
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip17]
  #   type = ElementAverageValue
  #   variable = slip17
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip18]
  #   type = ElementAverageValue
  #   variable = slip18
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip19]
  #   type = ElementAverageValue
  #   variable = slip19
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip20]
  #   type = ElementAverageValue
  #   variable = slip20
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip21]
  #   type = ElementAverageValue
  #   variable = slip21
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip22]
  #   type = ElementAverageValue
  #   variable = slip22
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip23]
  #   type = ElementAverageValue
  #   variable = slip23
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip24]
  #   type = ElementAverageValue
  #   variable = slip24
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip25]
  #   type = ElementAverageValue
  #   variable = slip25
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip26]
  #   type = ElementAverageValue
  #   variable = slip26
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip27]
  #   type = ElementAverageValue
  #   variable = slip27
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip28]
  #   type = ElementAverageValue
  #   variable = slip28
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip29]
  #   type = ElementAverageValue
  #   variable = slip29
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip30]
  #   type = ElementAverageValue
  #   variable = slip30
  #   #point = '1. 1. 1.'
  # [../]
  # [./slip31]
  #   type = ElementAverageValue
  #   variable = slip31
  #   #point = '1. 1. 1.'
  # [../]
  [./stress_center_xx]
    type = ElementAverageValue
    variable = stress_xx
    # #point = '1. 1. 1.'
    # use_displaced_mesh = true
  [../]
  [./stress_center_yy]
    type = ElementAverageValue
    variable = stress_yy
    #point = '1. 1. 1.'
  [../]
  [./stress_center_zz]
    type = ElementAverageValue
    variable = stress_zz
    #point = '1. 1. 1.'
  [../]
  [./strain_center_xx]
    type = ElementAverageValue
    variable = strain_xx
    #point = '1. 1. 1.'
  [../]
  [./stress_center_xy]
    type = ElementAverageValue
    variable = stress_xy
    #point = '1. 1. 1.'
  [../]
  [./strain_center_xy]
    type = ElementAverageValue
    variable = strain_xy
    #point = '1. 1. 1.'
  [../]
  [./u_xx]
    type = ElementAverageValue
    variable = disp_x
    #point = '1. 1. 1.'
  [../]
  [./u_yy]
    type = ElementAverageValue
    variable = disp_y
    #point = '1. 1. 1.'
  [../]
  [./u_zz]
    type = ElementAverageValue
    variable = disp_z
    #point = '1. 1. 1.'
  [../]
  [./pk2_xx]
    type = ElementAverageValue
    variable = pk2_xx
    #point = '1. 1. 1.'
  [../]
  [./pk2_yy]
    type = ElementAverageValue
    variable = pk2_yy
    #point = '1. 1. 1.'
  [../]
  [./pk2_zz]
    type = ElementAverageValue
    variable = pk2_zz
    #point = '1. 1. 1.'
  [../]
[]
