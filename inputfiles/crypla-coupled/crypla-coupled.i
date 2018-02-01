
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 50
  nz = 1
  xmin = 0
  xmax = 50 #[nm]
  ymin = -300
  ymax = 325
  zmin = 0
  zmax = 12.5
  elem_type = HEX8
  #uniform_refine = 3 # Initial uniform refinement of the mesh
  #parallel_type = REPLICATED
  #skip_partitioning = true # all nodes knows of the whole mesh
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  #enable_jit = true
  displacements = 'disp_x disp_y disp_z'
[]
[BCs]
  [./Periodic]
    [./x]
      auto_direction = 'x z'
      variable = 'c w c0 c1 c2 eta0 eta1 eta2'
    [../]
  [../]
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
  [./disp_z]
    type = PresetBC
    variable = disp_z
    boundary = 'back front'
    value = 0
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
    #initial_condition = 0
    #scaling = 1e3
  [../]

  # phase concentration  Sn in Cu
  [./c0]
    order = FIRST
    family = LAGRANGE
    #initial_condition = 0.0234
    initial_condition = 0.02
    #initial_condition = 0.1684 #makes df_cu/dx_cu=0
    #initial_condition = 0.10569
    #scaling = 1e2
  [../]
  # phase concentration  Sn in Cu3Sn
  [./c1]
    order = FIRST
    family = LAGRANGE
    #initial_condition = 0.2484
    #initial_condition = 0.2433
    initial_condition = 0.4350 #makes df_eps/dx_eps=0
    #scaling = 1e2
  [../]
  # phase concentration  Sn in Cu6Sn5
  [./c2]
    order = FIRST
    family = LAGRANGE
    #initial_condition = 0.435
    #initial_condition = 0.4351
    initial_condition = 0.9796 #makes df_eta/dx_eta=0
    #scaling = 1e2
  [../]
  # phase concentration  Sn in Sn

  # order parameter Cu
  [./eta0]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e3
  [../]
  # order parameter Cu6Sn5
  [./eta1]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
    #scaling = 1e3
  [../]
  # order parameter Sn
  [./eta2]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e3
  [../]


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
[ICs]
  [./eta_cu] #Cu
    type = FunctionIC
    variable = eta0
    #function = 'if(y<=0,1,0)'
    #function = '1-0.5*(1+tanh(y/100))' #close to equlibrium shape
    function = '1-0.5*(1+tanh(y/25))' #close to equlibrium shape
  [../]
  [./eta_sn] #Sn
    type = UnitySubVarIC
    variable = eta2
    etas = eta0
  [../]

  [./c] #Concentration of Sn
    type = VarDepIC
    variable = c
    cis = 'c0 c1 c2'
    etas = 'eta0 eta1 eta2'
  [../]
  #[./w]
  #  type = FunctionIC
  #  variable = w
  #  function = '-2635358969.*(1-0.5*(1+tanh(y/20)))+213062709.*(1-(1-0.5*(1+tanh(y/20))))'
  #[../]
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
    #prop_values = '0.5 200e-9 5e-10 1.5 16.29e-6 0' #J/m^2 m - ?
    prop_values = '0.5 50e-9 5e-10 1.5 16.29e-6 0' #J/m^2 m - ?
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
    outputs = exodus
    output_properties = mu
  [../]

  #Free energy
  [./fch_cu] #Chemical energy Cu phase
      type = DerivativeParsedMaterial
      f_name = fch0
      args = 'c0'
      material_property_names = 'A_cu B_cu C_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c0-chat_cu)^2+B_cu*(c0-chat_cu)+C_cu)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
  [../]
  [./fch_imc] #Chemical energy Cu6Sn5 phase grain 2
      type = DerivativeParsedMaterial
      f_name = fch1
      args = 'c1'
      material_property_names = 'A_eta B_eta C_eta chat_eta length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eta*(c1-chat_eta)^2+B_eta*(c1-chat_eta)+C_eta)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
  [../]
  [./fch_sn] #Chemical energy Sn phase
      type = DerivativeParsedMaterial
      f_name = fch2
      args = 'c2'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c2-chat_sn)^2+B_sn*(c2-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
  [../]
  #SwitchingFunction
  [./h_cu]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta0
      outputs = exodus
      output_properties = h0
      use_displaced_mesh = true
  [../]

  [./h_imc1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta1
      outputs = exodus
      output_properties = h1
      use_displaced_mesh = true
  [../]
  [./h_sn]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2'
      phase_etas = eta2
      outputs = exodus
      output_properties = h2
      use_displaced_mesh = true
  [../]

  #Double well
  [./g0]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta0
    well_only = True
    function_name = g0
    use_displaced_mesh = true
  [../]
  #Double well, not used
  [./g1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta1
    well_only = True
    function_name = g1
    use_displaced_mesh = true
  [../]

  [./g2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta2
    well_only = True
    function_name = g2
    use_displaced_mesh = true
  [../]

  [./Mgb]
    type=ParsedMaterial
    material_property_names = 'D_gb delta delta_real h0(eta0,eta1,eta2) h1(eta0,eta1,eta2) h2(eta0,eta1,eta2) A_cu A_eta A_sn length_scale energy_scale time_scale'
    f_name = Mgb
    function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h1*A_eta+h2*A_sn)*delta)'
    #function = '4e-5'
    use_displaced_mesh = true
  [../]
  [./CHMobility]
      type = DerivativeParsedMaterial
      f_name = M
      args = 'eta0 eta1 eta2'
      material_property_names = 'h0(eta0,eta1,eta2) h1(eta0,eta1,eta2) h2(eta0,eta1,eta2) D_cu D_eta D_sn A_cu A_eta A_sn Mgb length_scale energy_scale time_scale'
      #function = 's:=eta_cu^2+eta_imc1^2+eta_imc2^2+eta_sn^2;p:=eta_imc1^2*eta_imc2^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc1*h_imc2>1./16.,0,Mgb)' #nm^5/eVs
      function = '(length_scale^5/(energy_scale*time_scale))*(h0*D_cu/A_cu+h1*D_eta/A_eta+h2*D_sn/A_sn)' #'+h_imc1*h_imc2*Mgb' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
      derivative_order = 2
      #outputs = exodus
      use_displaced_mesh = true
  [../]

  [./ACMobility]
      type = GenericConstantMaterial
      prop_names = L
      prop_values = 2.7
  [../]

  [./time]
    type = TimeStepMaterial
    prop_time = time
    prop_dt = dt
    use_displaced_mesh = true
  [../]
  [./noise_constants]
    type = GenericConstantMaterial
    prop_names = 'T kb lambda dim' #temperature Boltzmann gridsize dimensionality
    prop_values = '493 8.6173303e-5 12.5 3'
  [../]
  [./nuc_eta]
    type =  DerivativeParsedMaterial
    f_name = nuc_eta
    #args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'time dt T kb lambda dim L h0 h2' #'h0(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3)'
    #function = '20*h0*h2'
    #function = 'sqrt(2*kb*T*L/(lambda^dim*dt))*h0*h3' #expression from Shen (2007) (without h0*h2)
    function = 'if(h0*h2>0.09&time<1,sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007) (without h0*h2)
    #function = '100'
    outputs = exodus
    use_displaced_mesh = true
  [../]

  # Elasticity
  [./elasticity_tensor_eta]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 701. #112.3 GPa in eV/nm^3
    poissons_ratio = 0.31
    base_name = eta1
  [../]
  [./strain_eta]
    type = ComputeFiniteStrain
    base_name = eta1
    eigenstrain_names = eT_eta
  [../]
  [./stress_eta]
    type = ComputeFiniteStrainElasticStress
    base_name = eta1
  [../]
  [./eigenstrain_eta]
    type = ComputeVariableEigenstrain
    #args = 'y'
    args = ' '
    base_name = eta1
    eigen_base = '1 1 1 0 0 0'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eta
    prefactor = pre
  [../]
  [./pre]
    type = DerivativeParsedMaterial
    material_property_names = 'h1'
    function = '0.01*h1'
    f_name = pre
    outputs = exodus
  [../]
  #[./pre]
  #  type = GenericConstantMaterial
  #  prop_names = pre
  #  #prop_values = 0.02
  #  prop_values = 0.01
  #  #prop_values = -0.003
  #[../]
  [./fel_eta]
    type = ElasticEnergyMaterial #This should use pk2 for large deformation
    args = ' '
    base_name = eta1
    f_name = fel1
    outputs = exodus
    output_properties = fel1
    use_displaced_mesh = true
  [../]

  [./elasticity_tensor_cu]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 936. #150 GPa in eV/nm^3
    poissons_ratio = 0.35
    base_name = eta0
  [../]
  [./strain_cu]
    type = ComputeFiniteStrain
    base_name = eta0
  [../]
  [./stress_cu]
    type = ComputeFiniteStrainElasticStress
    base_name = eta0
  [../]
  [./fel_cu]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta0
    f_name = fel0
    outputs = exodus
    output_properties = fel0
    use_displaced_mesh = true
  [../]

  [./crysp]
    type = FiniteStrainUObasedCPBaseName
    block = 0
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    base_name = 'eta2'
    maximum_substep_iteration = 5
  [../]
  [./strain2]
    type = ComputeFiniteStrain
    block = 0
    base_name = 'eta2'
  [../]
  [./elasticity_tensor2]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    block = 0
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 149.8022 137.31 137.31' #eV/nm^3 #From Darbandi 2013 table III
    fill_method = symmetric9
    euler_angle_1 = 120
    euler_angle_2 = 90
    euler_angle_3 = 0
    base_name = 'eta2'
  [../]
  #[./stress_sn]
  #  type = ComputeFiniteStrainElasticStress
  #  base_name = eta3
  #[../]
  [./fel_sn]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta2
    f_name = fel2
    outputs = exodus
    output_properties = fel2
    use_displaced_mesh = true
  [../]

  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'eta0 eta1 eta2'
    h          = 'h0 h1 h2'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeFiniteStrain
    base_name = global
  [../]

  # Total elastic energy
  [./Total_elastic_energy]
    type = DerivativeMultiPhaseMaterial
    f_name = fel_tot
    W = 0
    fi_names = 'fel0 fel1 fel2'
    hi_names = 'h0 h1 h2'
    etas = 'eta0 eta1 eta2'
    outputs = exodus
    output_properties = fel_tot
    g = g0
    use_displaced_mesh = true
  [../]

  #sum chemical and elastic energies
  [./F_cu]
    type = DerivativeSumMaterial
    f_name = F0
    args = 'c0'
    sum_materials = 'fch0 fel0'
    #sum_materials = 'fch0'
    use_displaced_mesh = true
  [../]
  [./F_eta]
    type = DerivativeSumMaterial
    f_name = F1
    args = 'c1'
    sum_materials = 'fch1 fel1'
    #sum_materials = 'fch1'
    use_displaced_mesh = true
  [../]
  [./F_sn]
    type = DerivativeSumMaterial
    f_name = F2
    args = 'c2'
    sum_materials = 'fch2 fel2'
    #sum_materials = 'fch2'
    use_displaced_mesh = true
  [../]

[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
    base_name = global
    use_displaced_mesh = true
    strain = FINITE
  [../]
  #Kernels for split Cahn-Hilliard equation
  # Cahn-Hilliard Equation
  [./CHBulk] # Gives the residual for the concentration, dF/dc-mu
      type = KKSSplitCHCRes
      variable = c
      ca       = c1
      cb       = c2
      fa_name  = F1 #only fa is used
      fb_name  = F2
      w        = w
      h_name   = h1
      args_a = ' '
      use_displaced_mesh = true
  [../]

  [./dcdt] # Gives dc/dt
      type = CoupledTimeDerivative
      variable = w
      v = c
      use_displaced_mesh = true
  [../]
  [./ckernel] # Gives residual for chemical potential dc/dt+M\grad(mu)
      type = SplitCHWRes
      mob_name = M
      variable = w
      args = 'eta0 eta1 eta2'
      use_displaced_mesh = true
  [../]

  #KKS conditions
  # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
  [./chempot_cu_eta]
    type = KKSPhaseChemicalPotential
    variable = c0
    cb       = c1
    fa_name  = F0
    fb_name  = F1
    use_displaced_mesh = true
  [../]
  [./chempot_eta_sn]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb       = c2
    fa_name  = F1
    fb_name  = F2
    use_displaced_mesh = true
  [../]
  [./phaseconcentration] # enforce c = sum h_i*c_i
    type = KKSMultiPhaseConcentration
    variable = c2
    cj = 'c0 c1 c2'
    hj_names = 'h0 h1 h2'
    etas = 'eta0 eta1 eta2'
    c = c
    use_displaced_mesh = true
  [../]

  #Kernels for Allen-Cahn equations
  [./KKSMultiACKernel]
    op_num = 3
    op_name_base = 'eta'
    ci_name_base = 'c'
    f_name_base = F
    #wi = 0.0624
    wi = 0.
    #wi = 4.
    #wi = 1.
    g_name_base = g
    use_displaced_mesh = true
  [../]

  #Nucleation Kernel
  #[./nucleation_eps]
  #  type = LangevinNoise # TODO: This draws random number from a uniform distribution, it should be from a standard gaussian instead(?) maybe always >0
  #  variable = eta1
  #  amplitude = 1
  #  seed = 123456789
  #  multiplier = nuc_eps
  #[../]
  #
  ##Nucleation Kernel
  #[./nucleation_eta]
  #  type = LangevinNoise
  #  variable = eta2
  #  amplitude = 1
  #  seed = 987654321
  #  multiplier = nuc_eta
  #[../]
  #Nucleation Kernel
  [./nucleation_eta]
    type = LangevinNoisePositive # TODO: This draws random number from a uniform distribution, it should be from a standard gaussian instead(?) maybe always >0
    variable = eta1
    amplitude = 1
    seed = 123456789
    multiplier = nuc_eta
    use_displaced_mesh = true
  [../]

[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
    flowprops = '1 32 0.001 0.05' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss
    base_name = 'eta2'
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    groups = '0 32'
    group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    #group_values = '23000' # MPa to make everything elastic
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 32
    #hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    hprops = '1.4 0.624 0.250 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]


[AuxVariables]
  # Visualization
  [./phase_id]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_density] #local free energy density
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_int]
      order = CONSTANT
      family = MONOMIAL
  [../]
  [./von_mises]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./von_mises_cu]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./von_mises_eta]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./von_mises_sn]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]

  [./von_mises_mpa]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./von_mises_cu_mpa]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./von_mises_eta_mpa]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./von_mises_sn_mpa]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_stress]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./hyd_stress_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_stress_cu]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./hyd_stress_cu_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_stress_eta]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./hyd_stress_eta_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_stress_sn]
    order = CONSTANT
    family = MONOMIAL
    outputs = none
  [../]
  [./hyd_stress_sn_mpa]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./hyd_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]

  [./phase_id]
    type = ParsedAux
    variable = phase_id
    args = 'h0 h1 h2'
    function = 'if(h0>0.5,1,0)+if(h1>0.5,2,0)+if(h2>0.5,3,0)'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]

  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor =global_stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./von_mises_kernel_cu]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises_cu
    rank_two_tensor =eta0_stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./von_mises_kernel_eta]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises_eta
    rank_two_tensor =eta1_stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./von_mises_kernel_sn]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises_sn
    rank_two_tensor =eta2_stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]

  [./hydrostatic_strain]
    #Calculates the global hydrostatic stress
    type = RankTwoScalarAux
    variable = hyd_strain
    rank_two_tensor =global_total_strain
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]
  [./hydrostatic_stress]
    #Calculates the global hydrostatic stress
    type = RankTwoScalarAux
    variable = hyd_stress
    rank_two_tensor =global_stress
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]
  [./hydrostatic_stress_cu]
    #Calculates the hydrostatic stress
    type = RankTwoScalarAux
    variable = hyd_stress_cu
    rank_two_tensor =eta0_stress
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]
  [./hydrostatic_stress_eta]
    #Calculates the hydrostatic stress
    type = RankTwoScalarAux
    variable = hyd_stress_eta
    rank_two_tensor =eta1_stress
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]
  [./hydrostatic_stress_sn]
    #Calculates the hydrostatic stress
    type = RankTwoScalarAux
    variable = hyd_stress_sn
    rank_two_tensor =eta2_stress
    execute_on = timestep_end
    scalar_type = Hydrostatic
  [../]

  # Convert stresses from eV/nm^3 to MPA
  [./hyd_mpa]
    type = ParsedAux
    variable = hyd_stress_mpa
    args = hyd_stress
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'hyd_stress*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./hyd_cu_mpa]
    type = ParsedAux
    variable = hyd_stress_cu_mpa
    args = 'hyd_stress_cu h0'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h0*hyd_stress_cu*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./hyd_eta_mpa]
    type = ParsedAux
    variable = hyd_stress_eta_mpa
    args = 'hyd_stress_eta h1'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h1*hyd_stress_eta*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./hyd_sn_mpa]
    type = ParsedAux
    variable = hyd_stress_sn_mpa
    args = 'hyd_stress_sn h2'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h2*hyd_stress_sn*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./vm_mpa]
    type = ParsedAux
    variable = von_mises_mpa
    args = von_mises
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'von_mises*to_MPa'
    execute_on = timestep_end
  [../]
  [./vm_cu_mpa]
    type = ParsedAux
    variable = von_mises_cu_mpa
    args = 'von_mises_cu h0'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h0*von_mises_cu*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./vm_eta_mpa]
    type = ParsedAux
    variable = von_mises_eta_mpa
    args = 'von_mises_eta h1'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h1*von_mises_eta*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./vm_sn_mpa]
    type = ParsedAux
    variable = von_mises_sn_mpa
    args = 'von_mises_sn h2'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h2*von_mises_sn*to_MPa'
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]



  [./f_density]
    type = KKSMultiFreeEnergy
    variable = f_density
    hj_names = 'h0 h1 h2'
    Fj_names = 'F0 F1 F2'
    gj_names = 'g0 g1 g2'
    additional_free_energy = f_int
    interfacial_vars = 'eta0 eta1 eta2'
    kappa_names = 'kappa kappa kappa'
    #w = 0.0624
    #w = 4.
    #w = 1.
    w = 0.
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./f_int]
    type = ParsedAux
    variable = f_int
    args = 'eta0 eta1 eta2'
    constant_names = 'sigma delta gamma length_scale energy_scale'
    constant_expressions = '0.5 60e-9 1.5 1e9 6.24150943e18'
    function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta0^4-0.5*eta0^2+0.25*eta1^4-0.5*eta1^2+0.25*eta2^4-0.5*eta2^2+gamma*(eta0^2*(eta1^2+eta2^2)+eta1^2*eta2^2)+0.25)'
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
[]
[Postprocessors]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  [./cu_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h0
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  [./imc1_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h1
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  [./sn_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h2
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  #Monitoring the progress
  [./step_size]
    type = TimestepSize
  [../]
[]
[Debug]
  show_var_residual_norms = true
  show_material_props = false
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  line_search = default
  #line_search = none
  #line_search = bt
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'
  l_max_its = 20
  nl_max_its = 20
  l_tol = 1.0e-4
  #l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-8 #1.0e-10
  nl_abs_tol = 1.0e-11#1.0e-11

  #num_steps = 2000
  end_time = 2000 #50 hours

  # Picard stuff
  #picard_abs_tol = 1e-8
  #picard_rel_tol = 1e-6
  #picard_max_its = 10

  #very simple adaptive time stepper
  scheme = implicit-euler
  #dtmax = 1

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
  #[./TimeStepper]
  #  type = SolutionTimeAdaptiveDT
  #  dt = 5e-3
  #[../]

  #[./Adaptivity]
  #  # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
  #  initial_adaptivity = 3 # Number of times mesh is adapted to initial condition
  #  refine_fraction = 0.7 # Fraction of high error that will be refined
  #  coarsen_fraction = 0.1 # Fraction of low error that will coarsened
  #  max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  #[../]
[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    solve_type = PJFNK
    full = true
    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          asm   nonzero'

    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -sub_pc_factor_shift_type -pc_factor_shift_type'
    petsc_options_value = '     121              1e-10     1e-8     1e-5          asm       ilu             nonzero                     nonzero'

    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type pc_factor_mat_solver_package'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu       ilu             nonzero             nonzero            superlu_dist'
    ##
    #petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
    #petsc_options_value = 'bjacobi       ilu            nonzero'

    #petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
    #petsc_options_value = 'asm       ilu            nonzero'

    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu    ilu   nonzero'
  [../]
  [./mydebug]
    type = FDP
    full = true
    petsc_options_iname = '-sub_pc_factor_shift_type -pc_factor_shift_type'
    petsc_options_value = 'nonzero nonzero'
  [../]

[]

[Outputs]
  file_base = 3d-crypla
  exodus = true
  #exodus = true
  csv = true
  print_perf_log = true
  print_linear_residuals = true
  interval = 1 #5
  print_mesh_changed_info = true

[]
