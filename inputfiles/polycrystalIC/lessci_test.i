FUNKAR INTE
[Mesh]
  type = GeneratedMesh
  dim = 2
  elem_type = QUAD4
  nx = 100
  ny = 150
  #nz = 1
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1500
  #zmin = 0
  #zmax = 10
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  #enable_jit = true
[]
[Variables]
  [./c]
    #scaling = 1e2
  [../]
  # chemical potential
  [./w]
    #scaling = 1e1
    initial_condition = 0
  [../]

  #Cu
  [./eta0]
    #scaling = 1e1
  [../]
  # Sn
  [./eta1]
    #initial_condition = 0
    #scaling = 1e1
  [../]
  [./eta2]
    #initial_condition = 1
    #scaling = 1e1
  [../]
  [./eta3]
    initial_condition = 0
    #scaling = 1e1
  [../]


  #Cu
  [./c0]
    initial_condition = 0.02
    #initial_condition = 0.1617
    #scaling = 1e2
  [../]
  #Sn
  [./c1]
    initial_condition = 0.9745
    #scaling = 1e2
  [../]
  #Cu6Sn5
  [./c3]
    initial_condition = 0.4350
    #scaling = 1e2
  [../]
[]
[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 3 # Number of grains
    variable = 'eta0 eta1 eta2'
    rand_seed = 25752
    #columnar_3D = true
  [../]
  #[./grain_tracker]
  #  type = GrainTracker
  #  threshold = 0.2
  #  connecting_threshold = 0.08
  #  compute_halo_maps = true # Only necessary for displaying HALOS
  #  variable = ' '
  #[../]
[]
[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      op_num = 3
      var_name_base = eta
      polycrystal_ic_uo = voronoi
    [../]
  [../]
  [./c] #Concentration of Sn
    type = VarDepIC
    variable = c
    cis = 'c0 c1 c1 c3'
    etas = 'eta0 eta1 eta2 eta3'
  [../]
  #[./w]
  #  type = FunctionIC
  #  variable = w
  #  #function = 'if(y<1000,-16.4486,-0.74034)' #abaqus
  #  function = 'if(y<0,-16.4486,-0.74034)' #generated
  #[../]
[]

[BCs]
  [./Periodic]
    #generated mesh
    [./xy]
      auto_direction = 'x'
      variable = 'eta0 eta1 eta2 eta3 c w c0 c1 c3'
    [../]
  [../]

[]

[Materials]
  [./time]
    type = TimeStepMaterial
    prop_time = time
    prop_dt = dt
    use_displaced_mesh = false
  [../]
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale time_scale'
    prop_values = '1e9 6.24150943e18 1.' #m to nm J to eV s to h
  [../]
  [./model_constants]
    type = GenericConstantMaterial
    prop_names = 'sigma delta delta_real gamma tgrad_corr_mult'
    prop_values = '0.5 40e-9 5e-10 1.5 0' #J/m^2 m - ?
    #prop_values = '0.5 60e-9 1.5 0' #J/m^2 m - ?
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
  [./h0]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta0
      use_displaced_mesh = false
      outputs = exodus
      output_properties = h0
  [../]
  [./h1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta1
      use_displaced_mesh = false
      outputs = exodus
      output_properties = h1
  [../]
  [./h2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta2
      use_displaced_mesh = false
      outputs = exodus
      output_properties = h2
  [../]
  [./h3]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta3
      use_displaced_mesh = false
      outputs = exodus
      output_properties = h3
  [../]
  [./g0]
    type = BarrierFunctionMaterial
    eta = eta0
    well_only = true
    function_name = g0
    g_order = SIMPLE
    use_displaced_mesh = false
  [../]
  [./g1]
    type = BarrierFunctionMaterial
    eta = eta1
    well_only = true
    function_name = g1
    g_order = SIMPLE
    use_displaced_mesh = false
  [../]
  [./g2]
    type = BarrierFunctionMaterial
    eta = eta2
    well_only = true
    function_name = g2
    g_order = SIMPLE
    use_displaced_mesh = false
  [../]
  [./g3]
    type = BarrierFunctionMaterial
    eta = eta3
    well_only = true
    function_name = g3
    g_order = SIMPLE
    use_displaced_mesh = false
  [../]
  [./ACMobility]
      type = GenericConstantMaterial
      prop_names = L
      prop_values = 2.7 #2.7
  [../]
  [./noise_constants]
    type = GenericConstantMaterial
    prop_names = 'T kb lambda dim' #temperature Boltzmann gridsize dimensionality
    #prop_values = '493 8.6173303e-5  3'
    prop_values = '493 8.6173303e-5 10 2'
  [../]
  [./nuc]
    type =  DerivativeParsedMaterial
    f_name = nuc
    material_property_names = 'time dt T kb lambda dim L h0(eta0,eta1,eta2,eta3) h1(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3)'
    function = 'if(time<1&h0*(h1+h2)>0.09,sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007)
    outputs = exodus
    output_properties = nuc
    use_displaced_mesh = false
  [../]
  #Constants The energy parameters are for 220 C
  [./energy_A]
    type = GenericConstantMaterial
    prop_names = 'A_cu A_eps A_eta A_sn'
    #prop_values = '1.0133e5/Vm 4e5/Vm 4.2059e6/Vm' #J/m^3
    prop_values = '1.7756e10 2.4555e11 2.4555e11 2.3033e10' #J/m^3 Aeps=Aeta=2e6
    #prop_values = '1.5929e10 2.4555e12 2.4555e12 2.3020e10' #J/m^3 Aeps = 2e7 Aeta = 2e7
  [../]
  [./energy_B]
    type = GenericConstantMaterial
    prop_names = 'B_cu B_eps B_eta B_sn'
    #prop_values = '-2.1146e4/Vm -6.9892e3/Vm 7.168e3/Vm' #J/m^3
    prop_values = '-2.6351e9 -1.4014e9 2.3251e7 2.14216e8' #J/m^3
    #prop_values = '-2.5789e9 -1.3733e9 2.3175e7 2.1406e8' #J/m^3
  [../]
  [./energy_C]
    type = GenericConstantMaterial
    prop_names = 'C_cu C_eps C_eta C_sn'
    #prop_values = '-1.2842e4/Vm -1.9185e4/Vm -1.5265e4/Vm' #J/m^3
    prop_values = '-1.1441e9 -1.7294e9 -1.7646e9 -1.646e9' #J/m^3
    #prop_values = '-1.1529e9 -1.7330e9 -1.7646e9 -1.646e9' #J/m^3
  [../]
  [./energy_c_ab]
    type = GenericConstantMaterial
    prop_names = 'c_cu_eps c_cu_eta c_cu_sn c_eps_cu c_eps_eta c_eps_sn c_eta_cu c_eta_eps c_eta_sn c_sn_cu c_sn_eps c_sn_eta'
    prop_values = '0.02 0.1957 0.6088 0.2383 0.2483 0.2495 0.4299 0.4343 0.4359 0.9789 0.9839 0.9889' #-
    #prop_values = '0.0234 0.198 0.6088 0.2479 0.2489 0.000 0.4345 0.4349 0.4351 0.9789 0.000 0.9889' #- Aeps = 2e7 Aeta = 2e7
  [../]
  [./energy_chat]
    type = GenericConstantMaterial
    prop_names = 'chat_cu chat_eps chat_eta chat_sn'
    prop_values = '0.02 0.2433 0.4351 0.9889' #-
    #prop_values = '0.0234 0.2484 0.4350 0.9889' #-
  [../]
  [./diffusion_constants]
    type = GenericConstantMaterial
    prop_names = 'D_cu D_eps D_eta D_sn'
    #prop_values = '1e-20 6e-16 1.5e-14 1e-13' # m^2/s #D12 best slightly slow
    #prop_values = '1e-20 9.5e-16 3e-14 1e-13' # m^2/s #D15
    prop_values = '1e-20 1.25e-15 3.1e-14 1e-13' # m^2/s #D16 BEST
    #prop_values = '1e-16 1.25e-15 3.1e-14 1e-13' # m^2/s #D17
    #outputs = exodus
  [../]
  [./D_gb]
    type = ParsedMaterial
    material_property_names = 'D_eta'
    f_name = D_gb
    function = '200*D_eta'
  [../]

  [./fch_cu] #Chemical energy cu phase
      type = DerivativeParsedMaterial
      f_name = fch0
      args = 'c0'
      material_property_names = 'A_cu B_cu C_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c0-chat_cu)^2+B_cu*(c0-chat_cu)+C_cu)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = false
      outputs = exodus
      output_properties = fch0
  [../]
  [./fch_sn] #Chemical energy Sn central grain
      type = DerivativeParsedMaterial
      f_name = fch1
      args = 'c1'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c1-chat_sn)^2+B_sn*(c1-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = false
      outputs = exodus
      output_properties = fch1
  [../]
  [./fch_imc] #Chemical energy Cu6Sn5
      type = DerivativeParsedMaterial
      f_name = fch3
      args = 'c3'
      material_property_names = 'A_eta B_eta C_eta chat_eta length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eta*(c3-chat_eta)^2+B_eta*(c3-chat_eta)+C_eta)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = false
      outputs = exodus
      output_properties = fch3
  [../]
  [./Mgb]
    type=ParsedMaterial
    material_property_names = 'D_gb delta delta_real h0(eta0,eta1,eta2,eta3) h1(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3) h3(eta0,eta1,eta2,eta3)  A_cu A_imc A_sn length_scale energy_scale time_scale'
    f_name = Mgb
    #function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h1*A_imc+(h2+h3)*A_sn)*delta)'
    function = 'if(h1*h2>0.09,(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h3*A_imc+(h2+h1)*A_sn)*delta),0)'
    #function = '4e-5'
    outputs = exodus
    output_properties = Mgb
  [../]
  [./CHMobility]
      type = DerivativeParsedMaterial
      f_name = M
      args = 'eta0 eta1 eta2 eta3'
      material_property_names = 'h0(eta0,eta1,eta2,eta3) h1(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3) h3(eta0,eta1,eta2,eta3) D_cu D_eta D_sn A_cu A_eta A_sn Mgb length_scale energy_scale time_scale'
      #function = 's:=eta_cu^2+eta_imc2^2+eta_imc1^2+eta_sn^2;p:=eta_imc2^2*eta_imc1^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc2*D_imc/A_imc+h_imc1*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc2*D_imc/A_imc+h_imc1*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc2*h_imc1>1./16.,0,Mgb)' #nm^5/eVs
      function = '(length_scale^5/(energy_scale*time_scale))*(h0*D_cu/A_cu+h3*D_eta/A_eta+(h2+h1)*D_sn/A_sn)+Mgb' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
      derivative_order = 2
      outputs = exodus
      output_properties = M
      use_displaced_mesh = false
  [../]

  [./F_cu]
    type = DerivativeSumMaterial
    f_name = F0
    args = 'c0'
    sum_materials = 'fch0'
    #sum_materials = 'fch2'
    use_displaced_mesh = false
  [../]
  [./F_eta]
    type = DerivativeSumMaterial
    f_name = F3
    args = 'c3'
    sum_materials = 'fch3'
    #sum_materials = 'fch1'
    use_displaced_mesh = false
  [../]
  [./F_sn1]
    type = DerivativeSumMaterial
    f_name = F1
    args = 'c1'
    sum_materials = 'fch1'
    #sum_materials = 'fch2'
    use_displaced_mesh = false
  [../]
  [./F_sn2]
    type = DerivativeSumMaterial
    f_name = F2
    args = 'c1'
    sum_materials = 'fch1'
    #sum_materials = 'fch2'
    use_displaced_mesh = false
  [../]
[]

[Kernels]
  #Nucleation of Cu6Sn5
  [./nuceta1]
    type = LangevinNoisePositive
    variable = eta3
    amplitude = 1
    seed = 12345
    multiplier = nuc
    use_displaced_mesh = false
  [../]
  # Cahn-Hilliard Equation
  [./CHBulk] # Gives the residual for the concentration, dF/dc-mu
      type = KKSSplitCHCRes
      variable = c
      ca       = c0
      cb       = c1
      fa_name  = F0 #only fa is used
      fb_name  = F1
      w        = w
      h_name   = h0
      args_a = 'eta0 eta1 eta2 eta3'
      use_displaced_mesh = false
  [../]

  [./dcdt] # Gives dc/dt
      type = CoupledTimeDerivative
      variable = w
      v = c
      use_displaced_mesh = false
  [../]
  [./ckernel] # Gives residual for chemical potential dc/dt+M\grad(mu)
      type = SplitCHWRes
      mob_name = M
      variable = w
      args = 'eta0 eta1 eta2 eta3'
      use_displaced_mesh = false
  [../]

  #KKS conditions
  # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
  [./chempot_cu_eta]
    type = KKSPhaseChemicalPotential
    variable = c0
    cb       = c1
    fa_name  = F0
    fb_name  = F1
    use_displaced_mesh = false
    args_a = ' '
    args_b = 'eta0 eta1 eta2 eta3'
  [../]
  [./chempot_eta_sn]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb       = c3
    fa_name  = F1
    fb_name  = F3
    use_displaced_mesh = false
    args_a = 'eta0 eta1 eta2 eta3'
  [../]
  [./phaseconcentration] # enforce c = sum h_i*c_i
    type = KKSMultiPhaseConcentration
    variable = c3
    cj = 'c0 c1 c1 c3'
    hj_names = 'h0 h1 h2 h3'
    etas = 'eta0 eta1 eta2 eta3'
    c = c
    use_displaced_mesh = false
  [../]

  [./KKSMultiACKernel]
    op_num = 4
    ci_num = 3
    op_name_base = 'eta'
    ci_name_base = 'c'
    f_name_base = F
    wi = 0.
    g_name_base = g
    ci_groups = '0 1 3 4'
    group_values = '0 1 3'
    use_displaced_mesh = false
  [../]
[]

[AuxVariables]
  [./f_density]
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
    hj_names = 'h0 h1 h2 h3'
    Fj_names = 'F0 F1 F2 F3'
    gj_names = 'g0 g1 g2 g3'
    additional_free_energy = f_int
    interfacial_vars = 'eta0 eta1 eta2 eta3'
    kappa_names = 'kappa kappa kappa kappa'
    #w = 0.0624
    #w = 4.
    #w = 1.
    w = 0.
    execute_on = 'initial timestep_end'
    use_displaced_mesh = false
  [../]
  [./f_int]
    type = ParsedAux
    variable = f_int
    args = 'eta0 eta1 eta2 eta3'
    constant_names = 'sigma delta gamma length_scale energy_scale'
    constant_expressions = '0.5 40e-9 1.5 1e9 6.24150943e18'
    #constant_expressions = '0.5 60e-9 1.5 1e9 6.24150943e18'
    function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta0^4-0.5*eta0^2+0.25*eta1^4-0.5*eta1^2+0.25*eta2^4-0.5*eta2^2+0.25*eta3^4-0.5*eta3^2+gamma*(eta0^2*(eta1^2+eta2^2+eta3^2)+eta1^2*(eta1^2+eta2^2+eta3^2)+eta2^2*(eta3^2))+0.25)'
    execute_on = 'initial timestep_end'
    use_displaced_mesh = false
  [../]
[]


[Postprocessors]
  [./imc_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h3
    execute_on = 'INITIAL TIMESTEP_END'
    use_displaced_mesh = false
  [../]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = false
  [../]
  [./step_size]
    type = TimestepSize
  [../]
[]
[Debug]
  show_var_residual_norms = true
  #show_material_props = true
[]
[Preconditioning]
  [./smp]
    type = SMP
    full = true #(Remove sone off-diagonal terms)
    solve_type = PJFNK

    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu       ilu             nonzero             nonzero            superlu_dist'
    #

    #petsc_options_iname = '-pc_asm_overlap -ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '1  121              1e-10     1e-8     1e-5          asm       ilu             nonzero             nonzero'
    #petsc_options_iname = '-pc_asm_overlap -ksp_gmres_restart  -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     1               121                asm       ilu             nonzero             nonzero'

    #Larrys suggestion
    petsc_options_iname = '-pc_asm_overlap -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type -sub_ksp_type'
    petsc_options_value = '2                  asm       lu             nonzero             nonzero                    preonly'

    #petsc_options_iname = '-pc_type -pc_factor_shift_type -sub_pc_factor_shift_type'
    #petsc_options_value = 'bjacobi       nonzero  nonzero'

    #petsc_options_iname = '-pc_type -pc_hypre_type   -pc_factor_shift_type  -ksp_gmres_restart'
    #petsc_options_value = 'hypre       boomeramg            nonzero    150'
  [../]


[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  #solve_type = 'NEWTON'

  #line_search = default
  #line_search = none
  #line_search = bt

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'

  l_max_its = 120
  nl_max_its = 25

  l_tol = 1.0e-3
  #l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-9 #1.0e-10
  nl_abs_tol = 1.0e-10#1.0e-11

  #num_steps = 2000
  end_time = 100 #50 hours
  #scheme = implicit-euler
  scheme = bdf2

  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 1e-3

      cutback_factor = 0.5
      growth_factor = 1.5
      optimal_iterations = 10
      linear_iteration_ratio = 25
      #postprocessor_dtlim = 5
  [../]

[]

[Outputs]
  file_base = 1Cu1imc2sn-lessci
  exodus = true
  csv = true
  print_perf_log = true
[]
