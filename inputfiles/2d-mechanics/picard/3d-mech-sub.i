[Mesh]
    type = GeneratedMesh
    dim = 3
    nx = 40
    ny = 60
    nz= 1
    xmin = 0
    xmax = 2000 #[nm]
    ymin = -1000
    ymax = 2000
    zmin = 0
    zmax = 50
    elem_type = HEX8
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  enable_jit = true
  displacements = 'disp_x disp_y disp_z'
[]
[Variables]
  # Displacements, solved by the sub app
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e3
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e3
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e3
  [../]
[]
[BCs]
  #[./Periodic]
  #  [./disp]
  #    auto_direction = 'x z'
  #    variable = 'disp_x disp_y disp_z'
  #  [../]
  #[./]
  [./disp_y]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom top'
    value = 0
  [../]
  [./disp_x]
    type = PresetBC
    variable = disp_x
    boundary = 'left right'
    value = 0
  [../]
  [./disp_z]
    type = PresetBC
    variable = disp_z
    boundary = 'front back'
    value = 0
  [../]
  #[./left_y]
  #  type = PresetBC
  #  variable = disp_y
  #  boundary = 'left'
  #  value = 0
  #[../]
  #[./left_x]
  #  type = PresetBC
  #  variable = disp_x
  #  boundary = 'left'
  #  value = 0
  #[../]
  #[./left_z]
  #  type = PresetBC
  #  variable = disp_z
  #  boundary = 'left'
  #  value = 0
  #[../]
[]
[AuxVariables]
  # order parameter Cu
  [./eta0]
      order = FIRST
      family = LAGRANGE
  [../]
  # order parameter for Cu6Sn5
  [./eta1]
      order = FIRST
      family = LAGRANGE
  [../]
  # order parameter Sn
  [./eta2]
      order = FIRST
      family = LAGRANGE
  [../]

[]

[Materials]
  #SwitchingFunction
  [./h_cu]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta0
  [../]

  [./h_imc1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta1
  [../]
  [./h_sn]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta2
  [../]

  [./elasticity_tensor_eta]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 701. #112.3 GPa in eV/nm^3
    poissons_ratio = 0.31
    base_name = eta1
  [../]
  [./strain_eta]
    type = ComputeSmallStrain
    base_name = eta1
    displacements = 'disp_x disp_y disp_z'
    eigenstrain_names = eT_eta
  [../]
  [./stress_eta]
    type = ComputeLinearElasticStress
    base_name = eta1
  [../]
  [./eigenstrain_eta]
    type = ComputeEigenstrain
    base_name = eta1
    eigen_base = '1 1 1'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eta
    prefactor = pre
  [../]
  [./pre]
    type = ParsedMaterial
    material_property_names = time
    f_name = pre
    #function = 'if(time<1,-0.003*time,-0.003)'
    function = '0.02'
  [../]


  [./elasticity_tensor_cu]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 936. #150 GPa in eV/nm^3
    poissons_ratio = 0.35
    base_name = eta0
  [../]
  [./strain_cu]
    type = ComputeSmallStrain
    base_name = eta0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_cu]
    type = ComputeLinearElasticStress
    base_name = eta0
  [../]

  [./elasticity_tensor_sn]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 119. #19 GPa in eV/nm^3
    poissons_ratio = 0.36
    base_name = eta2
  [../]
  [./strain_sn]
    type = ComputeSmallStrain
    base_name = eta2
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_sn]
    type = ComputeLinearElasticStress
    base_name = eta2
  [../]

  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'eta0 eta1 eta2'
    h          = 'h0 h1 h2'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y disp_z'
    base_name = global
  [../]

[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    base_name = global
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]
[Debug]
  show_var_residual_norms = true
  show_material_props = false
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  #solve_type = 'NEWTON'
  line_search = default
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  #petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  #petsc_options_value = 'asm lu 1 101'

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew '
  l_max_its = 30
  nl_max_its = 30
  l_tol = 1.0e-4
  l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-8 #1.0e-10
  nl_abs_tol = 1.0e-5#1.0e-11
  picard_rel_tol = 1e-8
  picard_abs_tol = 1e-9

  #num_steps = 2000
  end_time = 36000 #50 hours
  ##very simple adaptive time stepper
  #num_steps = 1
  scheme = bdf2
  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 0.01
      cutback_factor = 0.5
      growth_factor = 2.
      optimal_iterations = 15
  [../]
[]
[Outputs]
  exodus = true
  print_linear_residuals = false
[]
