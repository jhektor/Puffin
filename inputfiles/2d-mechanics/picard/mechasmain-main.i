[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 6
  ny = 25
  xmin = 0
  xmax = 600 #[nm]
  ymin = -1000
  ymax = 1500
  elem_type = QUAD4
  uniform_refine = 2 # Initial uniform refinement of the mesh
  parallel_type = REPLICATED
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  #enable_jit = true
  displacements = 'disp_x disp_y'
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
[]
[BCs]
  [./disp_y]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./disp_x]
    type = PresetBC
    variable = disp_x
    boundary = 'left right'
    value = 0
  [../]
[]
[AuxVariables]
  # order parameter Cu
  [./eta0]
      order = FIRST
      family = LAGRANGE
  [../]
  # order parameter for Cu3Sn
  [./eta1]
      order = FIRST
      family = LAGRANGE
  [../]
  # order parameter Cu6Sn5
  [./eta2]
      order = FIRST
      family = LAGRANGE
  [../]
  # sn
  [./eta3]
      order = FIRST
      family = LAGRANGE
  [../]

[]

[Materials]
  #SwitchingFunction
  [./h_cu]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta0
  [../]
  [./h_imc1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta1
  [../]
  [./h_imc2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta2
  [../]
  [./h_sn]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta3
  [../]

  # Elasticity
  [./elasticity_tensor_eta]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 701. #112.3 GPa in eV/nm^3
    poissons_ratio = 0.31
    base_name = eta2
  [../]
  [./strain_eta]
    type = ComputeSmallStrain
    base_name = eta2
    #displacements = 'disp_x disp_y'
    eigenstrain_names = eT_eta
  [../]
  [./stress_eta]
    type = ComputeLinearElasticStress
    base_name = eta2
  [../]
  [./eigenstrain_eta]
    type = ComputeEigenstrain
    base_name = eta2
    eigen_base = '1 1 0'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eta
    prefactor = pre
  [../]
  #[./pre]
  #  type = ParsedMaterial
  #  args = 'eta0 eta1 eta2 eta3'
  #  material_property_names = 'time h2(eta0,eta1,eta2,eta3)'
  #  f_name = pre
  #  #function = 'if(time<1,-0.003*time,-0.003)'
  #  #function = '-0.003'
  #  function = '0.02'
  #  #function = '0'
  #[../]
  [./pre]
    type = GenericConstantMaterial
    prop_names = pre
    prop_values = 0.02
  [../]
  [./fel_eta]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta2
    f_name = fel2
  [../]

  [./elasticity_tensor_eps]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 936. #150 GPa in eV/nm^3 WRONG VALUE
    poissons_ratio = 0.35
    base_name = eta1
  [../]
  [./strain_eps]
    type = ComputeSmallStrain
    base_name = eta1
    #displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_eps]
    type = ComputeLinearElasticStress
    base_name = eta1
  [../]
  [./fel_eps]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta1
    f_name = fel1
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
    #displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_cu]
    type = ComputeLinearElasticStress
    base_name = eta0
  [../]
  [./fel_cu]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta0
    f_name = fel0
  [../]

  [./elasticity_tensor_sn]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 119. #19 GPa in eV/nm^3
    poissons_ratio = 0.36
    base_name = eta3
  [../]
  [./strain_sn]
    type = ComputeSmallStrain
    base_name = eta3
    #displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress_sn]
    type = ComputeLinearElasticStress
    base_name = eta3
  [../]
  [./fel_sn]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta3
    f_name = fel3
  [../]

  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'eta0 eta1 eta2 eta3'
    h          = 'h0 h1 h2 h3'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeSmallStrain
    #displacements = 'disp_x disp_y disp_z'
  [../]

[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    base_name = global
    use_displaced_mesh = false
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  #solve_type = 'NEWTON'
  line_search = default
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'
  l_max_its = 100
  nl_max_its = 20
  l_tol = 1.0e-4
  #l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-6 #1.0e-10
  nl_abs_tol = 1.0e-8#1.0e-11

  #num_steps = 1
  #end_time = 180000 #50 hours

  #very simple adaptive time stepper
  scheme = implicit-euler
  #dtmax = 1
  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 5e-3

      cutback_factor = 0.5
      growth_factor = 1.25
      optimal_iterations = 2
      linear_iteration_ratio = 25
      #postprocessor_dtlim = 5
  [../]

  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    solve_type = PJFNK
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_factor_shift_type'
    petsc_options_value = '     121              1e-10     1e-8     1e-5          asm            nonzero'
    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu    ilu   nonzero'
    #petsc_options_iname = '-pc_type -pc_hypre_type -sub_pc_factor_shift_type'
    #petsc_options_value = 'hypre boomeramg nonzero'
  [../]
  [./mydebug]
    type = FDP
    full = true
  [../]
[]


[Outputs]
  file_base = mechmain
  [./exodus_out]
    type = Exodus
    interval =1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  print_linear_residuals = true
  interval = 1 #5
  print_mesh_changed_info = true
[]
[MultiApps]
  [./sub]
    type =  TransientMultiApp
    app_type = PuffinApp
    positions = '0 0 0'
    input_files = mechasmain-sub.i
    execute_on = TIMESTEP_BEGIN
    sub_cycling = true
    detect_steady_state = true
    steady_state_tol = 1e-08
  [../]
[]

[Transfers]
  [./dispx_to_sub]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_x
    variable = disp_x
  [../]
  [./dispy_from_sub]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_y
    variable = disp_y
  [../]
  [./eta0_to_sub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = eta0
    variable = eta0
  [../]
  [./eta1_to_sub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = eta1
    variable = eta1
  [../]
  [./eta2_to_sub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = eta2
    variable = eta2
  [../]
  [./eta3_to_sub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = eta3
    variable = eta3
  [../]
[]
