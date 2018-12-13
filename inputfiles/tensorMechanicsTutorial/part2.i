[Problem]
  coord_type = RZ #2D axisymmetric
[]
[GlobalParams]
  displacements = 'disp_r disp_z'
[]

[Mesh]
  file = necking_quad4.e
  uniform_refine = 0
  second_order = true #uses second order elements, understood by the master action
[]

[BCs]
  [./left]
    type = PresetBC
    variable = disp_r
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = PresetBC
    variable = disp_z
    boundary = bottom
    value = 0.0
  [../]
  [./top]
    type = FunctionPresetBC
    variable = disp_z
    boundary = top
    function = '0.0007*t'
  [../]
[]

#Tensor Mechanics action, sets up kernels and variables
[Modules/TensorMechanics/Master]
  [./block1]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy strain_yy stress_zz vonmises_stress'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeMultiPlasticityStress
    ep_plastic_tolerance = 1e-9
    plastic_models = J2
  [../]
[]

[UserObjects]
  [./str]
    type = TensorMechanicsHardeningCubic
    value_0 = 240
    value_residual = 300
    internal_0 = 0
    internal_limit = 0.005
  [../]
  [./J2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1e-3
    internal_constraint_tolerance = 1e-9
  [../]
[]

[Postprocessors]
  [./avg_stress_bottom]
    type = SideAverageValue
    variable = stress_yy
    boundary = bottom
  [../]
  [./avg_strain_bottom]
    type = SideAverageValue
    variable = strain_yy
    boundary = bottom
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  end_time = 40
  dt = 0.25
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
[]

[Outputs]
  exodus = true
  csv = true
  print_perf_log = true
  print_linear_residuals = false
[]
