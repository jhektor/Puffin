[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1000
  ny = 1
  xmin = 0
  xmax = 100000 # [nm]
[]

[BCs]
  [./neumann]
    type = NeumannBC
    boundary = 'left right'
    variable = 'c w c_cu c_imc c_sn eta_cu eta_imc eta_sn'
    value = 0
  [../]
  [./Periodic]
    [./y]
      auto_direction = y
      variable = 'c w c_cu c_imc c_sn eta_cu eta_imc eta_sn'
    [../]
  [../]
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
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c_cu]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  [../]
  [./c_imc]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.455
  [../]
  [./c_sn]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.95
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

[ICs]
  [./eta1]
    # Cu
    variable = eta_cu
    type = FunctionIC
    function = 'if(x<=50000,1,0)'
  [../]
  [./eta2]
    # Cu6Sn5
    # function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>8&r<=16,1,0)'
    variable = eta_imc
    type = FunctionIC
    function = 'if(x>50000&x<=58000,1,0)'
  [../]
  [./eta3]
    # Sn
    # function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>16,1,0)'
    variable = eta_sn
    type = FunctionIC
    function = if(x>58000,1,0)
  [../]
  [./c]
    # Concentration of Sn
    variable = c
    type = FunctionIC
    function = '0.05*if(x<=50000,1,0)+0.455*if(x>50000&x<=58000,1,0)+0.95*if(x>58000,1,0)'
  [../]
[]

[Materials]
  # scalings
  # Constants
  # Free energy
  # SwitchingFunction
  # Double well, not used
  # Double well, not used
  # Double well, not used
  # constant properties
  # [./CHMobility]
  # type = DerivativeParsedMaterial
  # f_name = M
  # args = 'eta_cu eta_imc eta_sn'
  # material_property_names = 'h_cu(eta_cu,eta_imc,eta_sn) h_imc(eta_cu,eta_imc,eta_sn) h_sn(eta_cu,eta_imc,eta_sn) D_cu D_imc D_sn A_cu A_imc A_sn length_scale energy_scale'
  # function = '(length_scale^5/energy_scale)*(h_cu*D_cu/A_cu+h_imc*D_imc/A_imc+h_sn*D_cu/A_sn)' #nm^5/eVs
  # [../]
  # [./ACMobility]
  # type = DerivativeParsedMaterial
  # f_name = L
  # args = 'eta_cu eta_imc eta_sn'
  # material_property_names = 'M(eta_cu,eta_imc,eta_sn)'
  # #function = 'M_cu*h_cu+M_imc*h_imc+M_sn*h_sn'
  # #function = '(0.5*eta_cu^2*eta_imc^2*(M_cu+M_imc)+0.5*eta_sn^2*eta_imc^2*(M_sn+M_imc)+0.5*eta_cu^2*eta_sn^2*(M_cu+M_sn))/(eta_cu^2*eta_imc^2+eta_cu^2*eta_sn^2+eta_sn^2*eta_imc^2)'
  # function = '0.5*M'
  # [../]
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale'
    prop_values = '1e9 6.24150943e18' # m to nm J to eV
  [../]
  [./constants]
    type = GenericConstantMaterial
    prop_names = 'A_cu A_imc A_sn D_cu D_imc D_sn sigma delta gamma tgrad_corr_mult'
    prop_values = '1e8 1e9 1e9 1e-25 1e-16 1e-14 0.5 0.667e-6 1.5 0' # J/m^3 m^2/s J/m^2 m - ?
  [../]
  [./kappa]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = kappa
    function = 0.75*sigma*delta*energy_scale/length_scale # eV/nm
  [../]
  [./mu]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = mu
    function = 6*sigma/delta*energy_scale/length_scale^3 # eV/nm^3
  [../]
  [./fch_cu]
    # Chemical energy Cu phase
    type = DerivativeParsedMaterial
    f_name = fch_cu
    args = c_cu
    material_property_names = 'A_cu length_scale energy_scale'
    function = (energy_scale/length_scale^3)*(0.5*A_cu*(c_cu-0.076)^2) # eV/nm^3
  [../]
  [./fch_imc]
    # Chemical energy Cu phase
    type = DerivativeParsedMaterial
    f_name = fch_imc
    args = c_imc
    material_property_names = 'A_imc length_scale energy_scale'
    function = (energy_scale/length_scale^3)*(0.5*A_imc*(c_imc-0.455)^2-1e6) # eV/nm^3
  [../]
  [./fch_sn]
    # Chemical energy Sn phase
    type = DerivativeParsedMaterial
    f_name = fch_sn
    args = c_sn
    material_property_names = 'A_sn length_scale energy_scale'
    function = (energy_scale/length_scale^3)*(0.5*A_sn*(c_sn-0.978)^2+1e7) # eV/nm^3
  [../]
  [./h_cu]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_cu
    all_etas = 'eta_cu eta_imc eta_sn'
    phase_etas = eta_cu
  [../]
  [./h_imc]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_imc
    all_etas = 'eta_cu eta_imc eta_sn'
    phase_etas = eta_imc
  [../]
  [./h_sn]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_sn
    all_etas = 'eta_cu eta_imc eta_sn'
    phase_etas = eta_sn
  [../]
  [./g_cu]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_cu
    well_only = True
    function_name = g_cu
  [../]
  [./g_imc]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_imc
    well_only = True
    function_name = g_imc
  [../]
  [./g_sn]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_sn
    well_only = True
    function_name = g_sn
  [../]
  [./mob]
    type = ParsedMaterial
    material_property_names = 'length_scale energy_scale D_imc A_imc'
    f_name = M
    function = (length_scale^5/energy_scale)*D_imc/A_imc
  [../]
  [./mobL]
    type = ParsedMaterial
    material_property_names = 'length_scale energy_scale D_imc A_imc D_sn A_sn sigma delta'
    f_name = L
    function = (length_scale^3/energy_scale)*(sigma/delta)*(D_sn/A_sn+D_imc/A_imc)/(3*0.75*sigma*delta*(0.999-0.476)^2)
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
    # Gives residual for chemical potential dc/dt+M\grad(mu)
    # args = 'eta_cu eta_imc eta_sn'
    type = SplitCHWRes # Ok if M is not depending on c or w
    mob_name = M
    variable = w
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
    # args      = 'eta_imc eta_sn'
    type = ACInterface
    variable = eta_cu
    kappa_name = kappa
    mob_name = L
    variable_L = false
  [../]
  [./ACdfintdeta_cu]
    # L*m*(eta_i^3-eta_i+2*beta*eta_i*sum_j eta_j^2)
    # args = 'eta_imc eta_sn'
    type = ACGrGrMulti # Jacobian not correct for non-constant mobility?
    variable = eta_cu
    v = 'eta_imc eta_sn'
    gamma_names = 'gamma gamma'
    mob_name = L
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
    # args      = 'eta_cu eta_sn'
    type = ACInterface
    variable = eta_imc
    kappa_name = kappa
    mob_name = L
    variable_L = false
  [../]
  [./ACdfintdeta_imc]
    # args = 'eta_cu eta_sn'
    type = ACGrGrMulti
    variable = eta_imc
    v = 'eta_cu eta_sn'
    gamma_names = 'gamma gamma'
    mob_name = L
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
    # args      = 'eta_cu eta_imc'
    type = ACInterface
    variable = eta_sn
    kappa_name = kappa
    mob_name = L
    variable_L = false
  [../]
  [./ACdfintdeta_sn]
    # args= 'eta_cu eta_imc'
    type = ACGrGrMulti
    variable = eta_sn
    v = 'eta_cu eta_imc'
    gamma_names = 'gamma gamma'
    mob_name = L
  [../]
[]

[AuxVariables]
  [./f_density]
    # local free energy density
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_int]
    order = CONSTANT
    family = MONOMIAL
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
    constant_names = 'gamma mu'
    constant_expressions = '1.5 1.'
    function = mu*(0.25*eta_cu^4-0.5*eta_cu^2+0.25*eta_imc^4-0.5*eta_imc^2+0.25*eta_sn^4-0.5*eta_sn^2+gamma*(eta_cu^2*(eta_imc^2+eta_sn^2)+eta_imc^2*eta_sn^2))+0.25
  [../]
[]

[Postprocessors]
  [./imc_thickness]
    type = NodalSum
    execute_on = timestep_end
    variable = eta_imc
  [../]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = TIMESTEP_END
  [../]
  [./imc_fraction]
    type = IMCFraction
    variable = eta_imc
    eta = 'eta_cu eta_sn'
    execute_on = Timestep_end
  [../]
[]

[Executioner]
  # very simple adaptive time stepper
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-11
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

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    full = true
  [../]
  [./mydebug]
    type = FDP
    full = true
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]

