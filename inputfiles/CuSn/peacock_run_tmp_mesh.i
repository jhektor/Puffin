# [Mesh]
# type = GeneratedMesh
# dim = 2
# nx = 2000
# ny = 1
# xmin = 0
# xmax = 200000 #[nm]
# ymin = 0
# ymax = 100
# elem_type = QUAD4
# []

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = 0
  xmax = 10000 # [nm]
  ymin = 0
  ymax = 10000
  elem_type = QUAD4
[]

[Variables]
  # concentration Sn
  # chemical potential
  # phase concentration  Sn in Cu
  # phase concentration  Sn in Cu6Sn5
  # phase concentration  Sn in Sn
  # order parameter Cu
  # order parameter Cu6Sn5
  # order parameter Sn
  [./c]
    # scaling = 1e-1
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    # scaling = 1e-3
    order = FIRST
    family = LAGRANGE
  [../]
  [./c_cu]
    # scaling = 1e-1
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  [../]
  [./c_imc]
    # scaling = 1e-2
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.455
  [../]
  [./c_sn]
    # scaling = 1e-5
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.95 # 0.999
  [../]
  [./eta_cu]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta_imc]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta_sn]
    # initial_condition = 0.0
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  # Kernels for split Cahn-Hilliard equation without composition gradent term(?)
  # Cahn-Hilliard Equation
  # 
  # KKS conditions
  # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
  # Kernels for Allen-Cahn equation for Cu
  # Kernels for Allen-Cahn equation for Cu6Sn5
  # Kernels for Allen-Cahn equation for Sn
  [./CHBulk]
    # Gives the residual for the concentration, dF/dc-mu
    # args_a = 'c_cu'
    type = KKSSplitCHCRes
    variable = c
    ca = c_imc
    cb = c_sn
    fa_name = fch_imc # only fa is used
    fb_name = fch_sn
    w = w
    h_name = h_imc
  [../]
  [./dcdt]
    # Gives dc/dt
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./ckernel]
    # Gives residual for chemical potential dc/dt+M\grad(mu) kansekse ska använda matDiffusion instaället, som kks_multiphase
    type = SplitCHWRes # Ok if M is not depending on c or w
    mob_name = M
    variable = w
    args = 'eta_cu eta_imc eta_sn'
  [../]
  [./chempot_cu_imc]
    type = KKSPhaseChemicalPotential
    variable = c_cu
    cb = c_imc
    fa_name = fch_cu
    fb_name = fch_imc
  [../]
  [./chempot_sn_cu]
    type = KKSPhaseChemicalPotential
    variable = c_imc
    cb = c_sn
    fa_name = fch_imc
    fb_name = fch_sn
  [../]
  [./phaseconcentration]
    # enforce c = sum h_i*c_i
    type = KKSMultiPhaseConcentration
    variable = c_sn
    cj = 'c_cu c_imc c_sn'
    hj_names = 'h_cu h_imc h_sn'
    etas = 'eta_cu eta_imc eta_sn'
    c = c
  [../]
  [./detadt_cu]
    type = TimeDerivative
    variable = eta_cu
  [../]
  [./ACBulkF_cu]
    # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used(?)
    type = KKSMultiACBulkF
    variable = eta_cu
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    gi_name = g_cu
    eta_i = eta_cu
    wi = 0
    mob_name = L
    args = 'c_cu c_imc c_sn eta_imc eta_sn'
  [../]
  [./ACBulkC_cu]
    # -L\sum_j dh_j/deta_i*mu_jc_j
    type = KKSMultiACBulkC
    variable = eta_cu
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    cj_names = 'c_cu c_imc c_sn'
    eta_i = eta_cu
    mob_name = L
    args = 'eta_imc eta_sn'
  [../]
  [./ACInterface_cu]
    # L*kappa*grad\eta_i
    type = ACInterface
    variable = eta_cu
    kappa_name = kappa
    mob_name = L
    args = 'eta_imc eta_sn'
    variable_L = true
  [../]
  [./ACdfintdeta_cu]
    # L*m*(eta_i^3-eta_i+2*beta*eta_i*sum_j eta_j^2)
    type = ACGrGrMulti # Jacobian not correct for non-constant mobility?
    variable = eta_cu
    v = 'eta_imc eta_sn'
    gamma_names = 'gamma gamma'
    mob_name = L
    args = 'eta_imc eta_sn'
  [../]
  [./detadt_imc]
    type = TimeDerivative
    variable = eta_imc
  [../]
  [./ACBulkF_imc]
    # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
    type = KKSMultiACBulkF
    variable = eta_imc
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    gi_name = g_imc
    eta_i = eta_imc
    wi = 0
    mob_name = L
    args = 'c_cu c_imc c_sn eta_cu eta_sn'
  [../]
  [./ACBulkC_imc]
    # -L\sum_j dh_j/deta_i*mu_jc_j
    type = KKSMultiACBulkC
    variable = eta_imc
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    cj_names = 'c_cu c_imc c_sn'
    eta_i = eta_imc
    mob_name = L
    args = 'eta_cu eta_sn'
  [../]
  [./ACInterface_imc]
    # L*kappa*grad\eta_i
    type = ACInterface
    variable = eta_imc
    kappa_name = kappa
    mob_name = L
    args = 'eta_cu eta_sn'
    variable_L = true
  [../]
  [./ACdfintdeta_imc]
    type = ACGrGrMulti
    variable = eta_imc
    v = 'eta_cu eta_sn'
    gamma_names = 'gamma gamma'
    mob_name = L
    args = 'eta_cu eta_sn'
  [../]
  [./detadt_sn]
    type = TimeDerivative
    variable = eta_sn
  [../]
  [./ACBulkF_sn]
    # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
    type = KKSMultiACBulkF
    variable = eta_sn
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    gi_name = g_sn
    eta_i = eta_sn
    wi = 0
    mob_name = L
    args = 'c_cu c_imc c_sn eta_imc eta_cu'
  [../]
  [./ACBulkC_sn]
    # -L\sum_j dh_j/deta_i*mu_jc_j
    type = KKSMultiACBulkC
    variable = eta_sn
    Fj_names = 'fch_cu fch_imc fch_sn'
    hj_names = 'h_cu h_imc h_sn'
    cj_names = 'c_cu c_imc c_sn'
    eta_i = eta_sn
    mob_name = L
    args = 'eta_cu eta_imc'
  [../]
  [./ACInterface_sn]
    # L*kappa*grad\eta_i
    type = ACInterface
    variable = eta_sn
    kappa_name = kappa
    mob_name = L
    args = 'eta_cu eta_imc'
    variable_L = true
  [../]
  [./ACdfintdeta_sn]
    type = ACGrGrMulti
    variable = eta_sn
    v = 'eta_cu eta_imc'
    gamma_names = 'gamma gamma'
    mob_name = L
    args = 'eta_cu eta_imc'
  [../]
[]

[AuxKernels]
  [./f_density]
    type = KKSMultiFreeEnergy
    variable = f_density
    hj_names = 'h_cu h_imc h_sn'
    Fj_names = 'fch_cu fch_imc fch_sn'
    gj_names = 'g_cu g_imc g_sn'
    additional_free_energy = f_int
    interfacial_vars = 'eta_cu eta_imc eta_sn'
    kappa_names = 'kappa kappa kappa'
    w = 0
  [../]
  [./f_int]
    type = ParsedAux
    variable = f_int
    args = 'eta_cu eta_imc eta_sn'
    constant_names = 'sigma delta gamma length_scale energy_scale'
    constant_expressions = '0.5 0.667e-6 1.5 1e9 6.24150943e18'
    function = 'mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta_cu^4-0.5*eta_cu^2+0.25*eta_imc^4-0.5*eta_imc^2+0.25*eta_sn^4-0.5*eta_sn^2+gamma*(eta_cu^2*(eta_imc^2+eta_sn^2)+eta_imc^2*eta_sn^2)+0.25)'
  [../]
[]

[Postprocessors]
  # Monitoring the progress
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = TIMESTEP_END
  [../]
  [./cu_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h_cu
    execute_on = TIMESTEP_END
  [../]
  [./imc_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h_imc
    execute_on = TIMESTEP_END
  [../]
  [./sn_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h_sn
    execute_on = TIMESTEP_END
  [../]
  [./time]
    type = RunTime
    time_type = active
  [../]
  [./step_size]
    type = TimestepSize
  [../]
[]

[Debug]
  show_var_residual_norms = false
  show_material_props = false
[]

[Executioner]
  # num_steps = 500
  # very simple adaptive time stepper
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-10 # 1.0e-11
  end_time = 64000000
  [./TimeStepper]
    # Turn on time stepping
    type = IterationAdaptiveDT
    dt = 1
    cutback_factor = 0.5
    growth_factor = 2.
    optimal_iterations = 5
  [../]
[]

[Outputs]
  # exodus = true
  file_base = keep/Cu-SnIC # moelans2011fig2_Limc_sn_Mvar_sharp
  csv = true
  print_perf_log = true
  interval = 1 # 5
  [./exodus_out]
    type = Exodus
    interval = 1 # 25
  [../]
[]

