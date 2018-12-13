[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  file = necking_quad4.e
  uniform_refine = 1
[]

[BCs]
  [./left]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./top]
    type = PresetBC
    variable = disp_y
    boundary = top
    value = 0.0035
  [../]
[]

#Tensor Mechanics action, sets up kernels and variables
[Modules/TensorMechanics/Master]
  [./block1]
    strain = SMALL #Seems to be plane strain by default
    add_variables = true
    generate_output = 'stress_xx strain_yy stress_zz strain_zz vonmises_stress'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options = -snes_ksp_ew
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
[]

[Outputs]
  exodus = true
  print_perf_log = true
[]
