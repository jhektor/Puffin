

[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 40
    ny = 60
    xmin = 0
    xmax = 2000 #[nm]
    ymin = 0
    ymax = 3000
    elem_type = QUAD4
[]

[BCs]
    [./neumann]
        type = NeumannBC
        boundary = 'top bottom'
        variable = 'c w c1 c2 c3 c0 eta1 eta2 eta3 eta0'
        value = 0
    [../]
    [./Periodic]
      [./x]
        auto_direction = 'x'
        variable = 'c w c1 c2 c3 c0 eta1 eta2 eta3 eta0'
      [../]
    [../]
[]

[Variables]
    # concentration Sn
    [./c]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e2
    [../]
    # chemical potential
    [./w]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]

    # phase concentration  Sn in Cu
    [./c1]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.002
        #initial_condition = 0.10569
        #scaling = 1e2
    [../]

    # phase concentration  Sn in Cu6Sn5
    [./c2]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.417
        #scaling = 1e2
    [../]
    # phase concentration  Sn in Cu6Sn5
    [./c3]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.417
        #scaling = 1e2
    [../]

    # phase concentration  Sn in Sn
    [./c0]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.999
        #scaling = 1e2
    [../]

    # order parameters
    [./PolycrystalVariables]
      op_num = 4
      var_name_base = eta
    [../]


[]
[ICs]
  [./eta_cu] #Cu
    type = UnitySubVarIC
    variable = eta1
    etas = 'eta2 eta3 eta0'
  [../]
  [./eta_imc1] #Cu6Sn5
    type = SmoothSuperellipsoidIC
    variable = eta2
    x1 = 1000
    y1 = 800
    n = 1.7
    a = 500
    b = 400
    int_width = 300
    invalue = 1
    outvalue = 0
  [../]
  [./eta_imc2] #Cu6Sn5
    type = SmoothSuperellipsoidIC
    variable = eta3
    x1 = 0
    y1 = 800
    n = 2
    a = 500
    b = 200
    int_width = 300
    invalue = 1
    outvalue = 0
  [../]
  [./eta_sn] #Sn
    type = UnitySubVarIC
    variable = eta0
    etas = 'eta2 eta3'
    use_threshold = true
    y_threshold = 800
  [../]

  [./c] #Concentration of Sn
    type = VarDepIC
    variable = c
    etas = 'eta1 eta2 eta3 eta0'
    cis = 'c1 c2 c3 c0'
  [../]
[]

[Materials]
  #scalings
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale time_scale'
    prop_values = '1e9 6.24150943e18 1.' #m to nm J to eV s to h
  [../]
  [./model_constants]
    type = GenericConstantMaterial
    prop_names = 'sigma delta delta_real gamma Vm tgrad_corr_mult'
    prop_values = '0.5 0.3e-6 5e-10 1.5 16.29e-6 0' #J/m^2 m - ?
  [../]
  #Constants
  [./energy_A]
    type = GenericConstantMaterial
    prop_names = 'A_cu A_imc A_sn'
    #prop_values = '1.0133e5/Vm 4e5/Vm 4.2059e6/Vm' #J/m^3
    prop_values = '6.2204e9 2.4555e10 2.5819e11' #J/m^3
  [../]
  [./energy_B]
    type = GenericConstantMaterial
    prop_names = 'B_cu B_imc B_sn'
    #prop_values = '-2.1146e4/Vm -6.9892e3/Vm 7.168e3/Vm' #J/m^3
    prop_values = '-1.2981e9 -0.0429e10 0.0044e11' #J/m^3
  [../]
  [./energy_C]
    type = GenericConstantMaterial
    prop_names = 'Cf_cu Cf_imc Cf_sn'
    #prop_values = '-1.2842e4/Vm -1.9185e4/Vm -1.5265e4/Vm' #J/m^3
    prop_values = '-0.7883e9 -0.1178e10 -0.0094e11' #J/m^3
    #prop_values = '0.9887e9 0. 8.37e8' #J/m^3
  [../]
  [./energy_c_ab]
    type = GenericConstantMaterial
    prop_names = 'c_imc_sn c_sn_imc'
    prop_values = '0.4529 0.9994' #-
  [../]
  [./energy_chat]
    type = GenericConstantMaterial
    prop_names = 'chat_cu chat_imc chat_sn'
    prop_values = '0.1057 0.41753 0.9994' #-
  [../]
  [./diffusion_constants]
    type = GenericConstantMaterial
    prop_names = 'D_cu D_imc D_sn'
    prop_values = '2.877e-36 6.575e-19 2.452e-17' # m^2/s
    #prop_values = '1e-19 6.575e-19 2.452e-17' # m^2/s
  [../]
  [./D_gb]
    type = ParsedMaterial
    material_property_names = 'D_imc'
    f_name = D_gb
    function = '200*D_imc'
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
  [../]

  [./L_imc0]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_imc A_sn A_imc c_imc_sn c_sn_imc length_scale energy_scale time_scale'
    f_name = L_imc0
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_imc/A_imc)/(3*kappa*(c_imc_sn-c_sn_imc)^2)'
  [../]

  #Free energy
  [./fch1] #Chemical energy Cu phase
      type = DerivativeParsedMaterial
      f_name = fch1
      args = 'c1'
      material_property_names = 'A_cu B_cu Cf_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c1-chat_cu)^2+B_cu*(c1-chat_cu)+Cf_cu)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch2] #Chemical energy Cu6Sn5 phase grain 1
      type = DerivativeParsedMaterial
      f_name = fch2
      args = 'c2'
      material_property_names = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c2-chat_imc)^2+B_imc*(c2-chat_imc)+Cf_imc)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch3] #Chemical energy Cu6Sn5 phase grain 2
      type = DerivativeParsedMaterial
      f_name = fch3
      args = 'c3'
      material_property_names = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c3-chat_imc)^2+B_imc*(c3-chat_imc)+Cf_imc)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch0] #Chemical energy Sn phase
      type = DerivativeParsedMaterial
      f_name = fch0
      args = 'c0'
      material_property_names = 'A_sn B_sn Cf_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c0-chat_sn)^2+B_sn*(c0-chat_sn)+Cf_sn)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  #SwitchingFunction
  [./h_cu]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta1 eta2 eta3 eta0'
      phase_etas = eta1
  [../]

  [./h_imc1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta1 eta2 eta3 eta0'
      phase_etas = eta2
  [../]
  [./h_imc2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta1 eta2 eta3 eta0'
      phase_etas = eta3
  [../]

  [./h_sn]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta1 eta2 eta3 eta0'
      phase_etas = eta0
  [../]

  #Double well, not used MAYBE USE TO KEEP THE ORDER PARAMETERS IN [0:1]
  [./g_cu]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta1
    well_only = True
    function_name = g1
  [../]
  #Double well, not used
  [./g_imc1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta2
    well_only = True
    function_name = g2
  [../]
  [./g_imc2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta3
    well_only = True
    function_name = g3
  [../]
  #Double well, not used
  [./g_sn]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta0
    well_only = True
    function_name = g0
  [../]

  [./CHMobility]
    type = DerivativeParsedMaterial
    f_name = M
    args = 'eta1 eta2 eta3 eta0'
    material_property_names = 'h1(eta1,eta2,eta3,eta0) h2(eta1,eta2,eta3,eta0) h3(eta1,eta2,eta3,eta0) h0(eta1,eta2,eta3,eta0) D_cu D_imc D_sn A_cu A_imc A_sn length_scale energy_scale time_scale Mgb'
    #function = 's:=eta1^2+eta2^2+eta3^2+eta0^2;p:=eta2^2*eta3^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc1*h_imc2>1./16.,0,Mgb)' #nm^5/eVs
    function = '(length_scale^5/(energy_scale*time_scale))*(h1*D_cu/A_cu+h2*D_imc/A_imc+h3*D_imc/A_imc+h0*D_sn/A_sn)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)' #'+h_imc1*h_imc2*(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h_cu*A_cu+h_imc1*A_imc+h_imc2*A_imc+h_sn*A_sn)*delta)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
    derivative_order = 2
    outputs = exodus_out
  [../]

  [./ACMobility]
    type = DerivativeParsedMaterial
    f_name = L
    args = 'eta1 eta2 eta3 eta0'
    material_property_names = 'L_cu_imc L_imc0 L_cu_sn L_imc_imc' # h_cu(eta1,eta_imc,eta0) h_imc(eta1,eta_imc,eta0) h_sn(eta1,eta_imc,eta0)'

    # Added epsilon to prevent division by 0 (Larry Aagesen)
    #function ='pf:=1e5;eps:=0.01;(L_cu_imc*(pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+L_imc0*((pf*eta2^2+eps)+(pf*eta3^2+eps))*(pf*eta0^2+eps)+L_cu_sn*(pf*eta1^2+eps)*(pf*eta0^2+eps)+L_imc_imc*(pf*eta2^2+eps)*(pf*eta3^2+eps))/((pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+((pf*eta2^2+eps)+(pf*eta3^2+eps))*(pf*eta0^2+eps)+(pf*eta1^2+eps)*(pf*eta0^2+eps)+(pf*eta2^2+eps)*(pf*eta3^2+eps))'
    #function ='pf:=1e5;eps:=1e-5;(L_cu_imc*(pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+L_imc0*((pf*eta2^2+eps)+(pf*eta3^2+eps))*(pf*eta0^2+eps)+L_cu_sn*(pf*eta1^2+eps)*(pf*eta0^2+eps)+L_imc_imc*(pf*eta2^2+eps)*(pf*eta3^2+eps))/((pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+((pf*eta2^2+eps)+(pf*eta3^2+eps))*(pf*eta0^2+eps)+(pf*eta1^2+eps)*(pf*eta0^2+eps)+(pf*eta2^2+eps)*(pf*eta3^2+eps))'
    function ='L_imc0'

    # Conditional function (Daniel Schwen)
    #function ='numer:=L_cu_imc*eta1^2*(eta2^2+eta3^2)+L_imc0*(eta2^2+eta3^2)*eta0^2+L_cu_sn*eta1^2*eta0^2;denom:=eta1^2*(eta2^2+eta3^2)+(eta2^2+eta3^2)*eta0^2+eta1^2*eta0^2;if(denom!=0,numer/denom,0.5*(L_cu_imc+L_imc0))'

    derivative_order = 2
    outputs = exodus_out
  [../]
[]

[Kernels]
    #Kernels for split Cahn-Hilliard equation
    # Cahn-Hilliard Equation
    [./CHBulk] # Gives the residual for the concentration, dF/dc-mu
        type = KKSSplitCHCRes
        variable = c
        ca       = c2
        cb       = c0
        fa_name  = fch2 #only fa is used
        fb_name  = fch0
        w        = w
        h_name   = h2
    [../]

    [./dcdt] # Gives dc/dt
        type = CoupledTimeDerivative
        variable = w
        v = c
    [../]
    [./ckernel] # Gives residual for chemical potential dc/dt+M\grad(mu)
        type = SplitCHWRes
        mob_name = M
        variable = w
        args = 'eta1 eta2 eta3 eta0'
    [../]

    #KKS conditions
    # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
    [./chempot_cu_imc]
      type = KKSPhaseChemicalPotential
      variable = c1
      cb       = c2
      fa_name  = fch1
      fb_name  = fch2
    [../]
    [./chempot_imc_imc]
      type = KKSPhaseChemicalPotential
      variable = c2
      cb       = c3
      fa_name  = fch2
      fb_name  = fch3
    [../]
    [./chempot_sn_cu]
      type = KKSPhaseChemicalPotential
      variable = c3
      cb       = c0
      fa_name  = fch3
      fb_name  = fch0
    [../]
    [./phaseconcentration] # enforce c = sum h_i*c_i
      type = KKSMultiPhaseConcentration
      variable = c0
      cj = 'c1 c2 c3 c0'
      hj_names = 'h1 h2 h3 h0'
      etas = 'eta1 eta2 eta3 eta0'
      c = c
    [../]

    #Kernels for Allen-Cahn equation for all order parameters
    [./KKSMultiACKernel]
      op_num = 4
      op_name_base = 'eta'
      ci_name_base ='c'
      wi = 10.
    [../]


[]

[AuxVariables]
    [./f_density] #local free energy density
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./f_int]
        order = CONSTANT
        family = MONOMIAL
    [../]
    [./s]
      order = FIRST
      family = LAGRANGE
    [../]
[]

[AuxKernels]
    [./f_density]
        type = KKSMultiFreeEnergy
        variable = f_density
        hj_names = 'h1 h2 h3 h0'
        Fj_names = 'fch1 fch2 fch3 fch0'
        gj_names = 'g1 g2 g3 g0'
        additional_free_energy = f_int
        interfacial_vars = 'eta1 eta2 eta3 eta0'
        kappa_names = 'kappa kappa kappa kappa'
        w = 10
        execute_on = 'initial timestep_end'
    [../]
    [./f_int]
        type = ParsedAux
        variable = f_int
        args = 'eta1 eta2 eta3 eta0'
        constant_names = 'sigma delta gamma length_scale energy_scale'
        constant_expressions = '0.5 0.4e-6 1.5 1e9 6.24150943e18'
        function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta1^4-0.5*eta1^2+0.25*eta2^4-0.5*eta2^2+0.25*eta3^4-0.5*eta3^2+0.25*eta0^4-0.5*eta0^2+gamma*(eta1^2*(eta2^2+eta3^2+eta0^2)+eta2^2*(eta3^2+eta0^2))+0.25)'
        execute_on = 'initial timestep_end'
    [../]
    [./s]
      type = ParsedAux
      variable = s
      args = 'eta1 eta2 eta3 eta0'
      #function = 'eta1^2*eta_imc^2+eta_imc^2*eta0^2+eta1^2*eta0^2'
      function = 'eta1+eta2+eta3+eta0'
    [../]
[]
[Postprocessors]
    [./total_energy]
      type = ElementIntegralVariablePostprocessor
      variable = f_density
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./cu_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h1
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc1_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h2
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc2_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h3
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./sn_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h0
      execute_on = 'Initial TIMESTEP_END'
    [../]
    #Monitoring the progress
    [./time]
      type = RunTime
      time_type = active
    [../]
    [./step_size]
      type = TimestepSize
    [../]
[]
[Debug]
  show_var_residual_norms = true
  show_material_props = false
  show_actions = false

  #show_parser = true

[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  line_search = default
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'


  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8 #1.0e-10
  nl_abs_tol = 1.0e-7#1.0e-11

  #num_steps = 2000
  end_time = 1e6
  #very simple adaptive time stepper
  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 1
      cutback_factor = 0.2
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
  file_base = kernelaction
  [./exodus_out]
    type = Exodus
    interval = 1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  interval = 1 #5
[]
