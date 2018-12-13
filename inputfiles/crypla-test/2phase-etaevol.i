# TODO: Elastic energy should be T*E^e not sigma*eps^e for finte strain? NOW IT'S T*E for Sn and sigma*eps for the rest
# TODO: Plastic energy is not included
# TODO: Check how eigenstrain is applied in finite strains: Eigenstrain are subtracted from strain increment, maybe should use F=F^eF^* as in plasticity?
# TODO: ComputeFiniteStrain does not compute Green-Lagrange strain.
[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 30
  ny = 30
  nz = 3
  xmin = -500
  xmax = 500
  ymin = -500
  ymax = 500
  zmin = 0
  zmax = 100
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./eta1]
    #initial_condition = 1
  [../]
  [./eta2]
    #initial_condition = 0
  [../]
[]

[ICs]
  [./eta1]
    type = SmoothCircleIC
    variable = eta1
    radius = 100
    x1 = 0
    y1 = 0
    3D_spheres = false
    invalue = 0
    outvalue = 1
    int_width = 120
  [../]
  [./eta2]
    type = UnitySubVarIC
    variable = eta2
    etas = eta1
  [../]
[]

[BCs]
  [./Periodic]
    [./x]
      auto_direction = 'x z'
      variable = 'eta1 eta2'
    [../]
  [../]
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
    group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    #group_values = '23000' # MPa to make everything elastic
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss1
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss1]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 32
    #hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    hprops = '1.4 0.624 0.250 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
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
    maximum_substep_iteration = 5
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
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 149.8022 137.31 137.31' #eV/nm^3 #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 60
    euler_angle_3 = 0
    base_name = 'eta1'
  [../]
  [./fe1]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe1
    base_name = eta1
    use_displaced_mesh = true #Not sure
    outputs = exodus
    output_properties = fe1
  [../]
  [./stress2]
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
    youngs_modulus = 700.92 #112.3 GPa
    poissons_ratio = 0.31
    base_name = eta2
  [../]
  [./eigenstrain_eps]
    type = ComputeVariableEigenstrain
    #args = ' '
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
    args = 'eta1 eta2'
    f_name = pre_eps
    material_property_names = 'h2(eta1,eta2) time'
    function = '0.02*h2*if(time>0.1,1,10*time)'
    outputs = exodus
    output_properties = pre_eps
  [../]
  [./fe2]
    type = ElasticEnergyMaterial
    args = 'eta1 eta2'
    f_name = fe2
    base_name = eta2
    use_displaced_mesh = true #Not sure
    outputs = exodus
    output_properties = fe2
  [../]
  #[./pre_eps]
  #  type = GenericFunctionMaterial
  #  prop_names = pre_eps
  #  prop_values = tdisp
  #  #prop_values = 0.02
  #  outputs = exodus
  #  output_properties = pre_eps
  #[../]
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale time_scale'
    prop_values = '1e9 6.24150943e18 1.' #m to nm J to eV s to h
  [../]
  [./model_constants]
    type = GenericConstantMaterial
    prop_names = 'sigma delta gamma tgrad_corr_mult'
    prop_values = '0.5 120e-9 1.5 0' #J/m^2 m - ?
  [../]
  [./kappa]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = kappa
    function = '0.75*sigma*delta*energy_scale/length_scale' #eV/nm
  [../]
  [./mu]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = mu
    function = '6*(sigma/delta)*energy_scale/length_scale^3' #eV/nm^3
    output_properties = mu
  [../]
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
  [./g1]
    type = BarrierFunctionMaterial
    eta = eta1
    well_only = true
    function_name = g1
    g_order = SIMPLE
    use_displaced_mesh = true
  [../]
  [./g2]
    type = BarrierFunctionMaterial
    eta = eta2
    well_only = true
    function_name = g2
    g_order = SIMPLE
    use_displaced_mesh = true
  [../]
  [./ACMobility]
      type = GenericConstantMaterial
      prop_names = L
      prop_values = 2.7
  [../]
  [./noise_constants]
    type = GenericConstantMaterial
    prop_names = 'T kb lambda dim' #temperature Boltzmann gridsize dimensionality
    prop_values = '493 8.6173303e-5 30 3'
  [../]
  [./nuc]
    type =  DerivativeParsedMaterial
    f_name = nuc
    material_property_names = 'time dt T kb lambda dim L' #'h0(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3)'
    function = 'if(time<1,sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007)
    outputs = exodus
    output_properties = nuc
    use_displaced_mesh = true
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

[Kernels]
  #Stress divergence
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    base_name = global
    strain = FINITE #this also sets incremental strain =true
  [../]

  #Nucleation of Cu6Sn5
  #[./nuceta2]
  #  type = LangevinNoisePositive
  #  variable = eta2
  #  amplitude = 1
  #  seed = 123456789
  #  multiplier = nuc
  #[../]

  # AC Kernels for eta1
  [./dt1]
    type = TimeDerivative
    variable = eta1
    use_displaced_mesh = true
  [../]
  # elastic energy contribution eta1
  [./dfedeta1]
    type = KKSMultiACBulkF
    variable = eta1
    args = 'eta2'
    eta_i = eta1
    Fj_names = 'fe1 fe2'
    hj_names = 'h1 h2'
    gi_name = 'g1'
    wi = 0
    mob_name = L
    use_displaced_mesh = true
  [../]
  # Gradient term eta1
  [./acint1]
    type = ACInterface
    variable = eta1
    mob_name = L
    kappa_name = kappa
    variable_L = false
    use_displaced_mesh = true
  [../]
  #Interface energy eta1
  [./finteta1]
    type = ACGrGrMulti
    variable = eta1
    mob_name = L
    gamma_names = 'gamma'
    v = eta2
    use_displaced_mesh = true
  [../]

  # AC Kernels for eta2
  [./dt2]
    type = TimeDerivative
    variable = eta2
    use_displaced_mesh = true
  [../]
  # elastic energy contribution eta2
  [./dfedeta2]
    type = KKSMultiACBulkF
    variable = eta2
    args = eta1
    eta_i = eta2
    Fj_names = 'fe1 fe2'
    hj_names = 'h1 h2'
    gi_name = 'g2'
    wi = 0
    mob_name = L
    use_displaced_mesh = true
  [../]
  # Gradient term eta2
  [./acint2]
    type = ACInterface
    variable = eta2
    mob_name = L
    kappa_name = kappa
    variable_L = false
    use_displaced_mesh = true
  [../]
  #Interface energy eta2
  [./finteta2]
    type = ACGrGrMulti
    variable = eta2
    mob_name = L
    gamma_names = 'gamma'
    v = eta1
    use_displaced_mesh = true
  [../]
[]

[AuxVariables]
  [./hyd_g]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./hyd_g_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_sn]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_sn_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_imc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_imc_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vm_g]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./vm_g_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #[./fe1]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  outputs = none
  #[../]
  #[./fe2]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  outputs = none
  #[../]
  [./hfe1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hfe2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  #[./fe1]
  #  type = MaterialRealAux
  #  variable = fe1
  #  property = fe1
  #  use_displaced_mesh = true
  #[../]
  #[./fe2]
  #  type = MaterialRealAux
  #  variable = fe2
  #  property = fe2
  #  use_displaced_mesh = true
  #[../]
  [./hfe1]
    type = ParsedAux
    variable = hfe1
    args = 'fe1 h1'
    function = 'h1*fe1'
    use_displaced_mesh = true
  [../]
  [./hfe2]
    type = ParsedAux
    variable = hfe2
    args = 'fe2 h2'
    function = 'h2*fe2'
    use_displaced_mesh = true
  [../]


  [./hyd_g]
    type = RankTwoScalarAux
    variable = hyd_g
    rank_two_tensor = global_stress
    scalar_type = Hydrostatic
    execute_on = timestep_end
  [../]
  [./hyd_g_mpa]
    type = ParsedAux
    variable = hyd_g_mpa
    args = hyd_g
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'hyd_g*to_MPa'
    execute_on = timestep_end
  [../]
  [./hyd_sn]
    type = RankTwoScalarAux
    variable = hyd_sn
    rank_two_tensor = eta1_stress
    scalar_type = Hydrostatic
    execute_on = timestep_end
  [../]
  [./hyd_sn_mpa]
    type = ParsedAux
    variable = hyd_sn_mpa
    args = 'h1 hyd_sn'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h1*hyd_sn*to_MPa'
    execute_on = timestep_end
  [../]
  [./hyd_imc]
    type = RankTwoScalarAux
    variable = hyd_imc
    rank_two_tensor = eta2_stress
    scalar_type = Hydrostatic
    execute_on = timestep_end
  [../]
  [./hyd_imc_mpa]
    type = ParsedAux
    variable = hyd_imc_mpa
    args = 'h2 hyd_imc'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h2*hyd_imc*to_MPa'
    execute_on = timestep_end
  [../]
  [./vm_g]
    type = RankTwoScalarAux
    variable = vm_g
    rank_two_tensor = global_stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./vm_g_mpa]
    type = ParsedAux
    variable = vm_g_mpa
    args = vm_g
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'vm_g*to_MPa'
    execute_on = timestep_end
  [../]
[]


[Postprocessors]
  [./imc_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h2
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  #[./e_yy1]
  #  type = PointValue
  #  variable = e_yy1
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./stress_yy1]
  #  type = PointValue
  #  variable = stress_yy1
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./e_yy2]
  #  type = PointValue
  #  variable = e_yy2
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./stress_yy2]
  #  type = PointValue
  #  variable = stress_yy2
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./e_yyg]
  #  type = PointValue
  #  variable = eyyg
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./stress_yyg]
  #  type = PointValue
  #  variable = syyg
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./stress_xyg]
  #  type = PointValue
  #  variable = stress_yy1
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]
  #[./stress_xzg]
  #  type = PointValue
  #  variable = sxzg
  #  point = '-250 0 50'
  #  #use_displaced_mesh = true
  #[../]

[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu       ilu             nonzero             nonzero            superlu_dist'

    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
    petsc_options_value = '     121              1e-10     1e-8     1e-5          asm       ilu             nonzero             nonzero            superlu_dist'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  line_search = default
  #line_search = none
  #line_search = bt
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'
  l_max_its = 100
  nl_max_its = 20
  l_tol = 1.0e-4
  #l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-8 #1.0e-10
  nl_abs_tol = 1.0e-10#1.0e-11

  #num_steps = 2000
  end_time = 2000 #50 hours
  scheme = implicit-euler

  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 5e-3

      cutback_factor = 0.5
      growth_factor = 1.25
      optimal_iterations = 10
      linear_iteration_ratio = 25
      #postprocessor_dtlim = 5
  [../]

[]

[Outputs]
  file_base = 2phase-eigenstrain-eVnm
  exodus = true
[]
