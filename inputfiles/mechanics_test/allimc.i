[Mesh]
    type = GeneratedMesh
    dim = 3
    nx = 10
    ny = 200
    nz = 1
    xmin = 0
    xmax = 500 #[nm]
    ymin = -5000
    ymax = 5000
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
#[Adaptivity]
#  marker = 'combo'
#
#  [./Indicators]
#    [./error1]
#      type = GradientJumpIndicator
#      variable = eta1
#    [../]
#    [./error2]
#      type = GradientJumpIndicator
#      variable = eta2
#
#    [../]
#  [../]
#
#  [./Markers]
#    [./combo]
#      type = ComboMarker
#      markers = 'errorfrac1 errorfrac2'
#    [../]
#    [./errorfrac1]
#      type = ErrorFractionMarker
#      indicator = error1
#      refine = 0.75
#      coarsen = 0.1
#    [../]
#    [./errorfrac2]
#      type = ErrorFractionMarker
#      indicator = error2
#      refine = 0.75
#      coarsen = 0.1
#    [../]
#  [../]
#[]

[BCs]
    [./Periodic]
      [./xz]
        auto_direction = 'x z'
        variable = 'c w c0 c1 c2 c3 eta0 eta1 eta2 eta3'
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
      boundary = 'bottom left right'
      value = 0
    [../]
    [./disp_z]
      type = PresetBC
      variable = disp_z
      boundary = 'bottom front back'
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
        initial_condition = 0.2490 #makes df_eps/dx_eps=0
        #scaling = 1e2
    [../]
    # phase concentration  Sn in Cu6Sn5
    [./c2]
        order = FIRST
        family = LAGRANGE
        #initial_condition = 0.435
        #initial_condition = 0.4351
        initial_condition = 0.4350 #makes df_eta/dx_eta=0
        #scaling = 1e2
    [../]

    # phase concentration  Sn in Sn
    [./c3]
        order = FIRST
        family = LAGRANGE
        #initial_condition = 0.95
        #initial_condition = 0.9889
        initial_condition = 0.9796 #makes df_sn/dx_sn=0
        #scaling = 1e2
    [../]

    # order parameter Cu
    [./eta0]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]
    # order parameter Cu3Sn
    [./eta1]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
        #scaling = 1e3
    [../]
    # order parameter Cu6Sn5
    [./eta2]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0
        #scaling = 1e3
    [../]

    # order parameter Sn
    [./eta3]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]

    [./disp_x]
      order = FIRST
      family = LAGRANGE
      #scaling = 1e-3
    [../]
    [./disp_y]
      order = FIRST
      family = LAGRANGE
      #scaling = 1e-3
    [../]
    [./disp_z]
      order = FIRST
      family = LAGRANGE
      #scaling = 1e-3
    [../]

[]
[ICs]
    [./eta_cu] #Cu
      type = FunctionIC
      variable = eta0
      #function = 'if(y<=0,1,0)'
      function = '1-0.5*(1+tanh(y/150))' #close to equlibrium shape

    [../]
    #[./eta_imc1] #Cu3Sn
    #  type = FunctionIC
    #  variable = eta1
    #  function = 'if(x>10000&x<=12600,1,0)'
    #[../]
    #[./eta_imc2] #Cu6Sn5
    #  type = FunctionIC
    #  variable = eta2
    #  function = 'if(x>10000&x<=19800,1,0)'
    #[../]
    [./eta_sn] #Sn
      type = UnitySubVarIC
      variable = eta3
      etas = eta0
      #function = 'if(x>19800,1,0)'
      #function = 'if(y>0,1,0)'
      #function = '0.5*(1+tanh(y/200))'
    [../]

    [./c] #Concentration of Sn
      type = VarDepIC
      variable = c
      cis = 'c0 c1 c2 c3'
      etas = 'eta0 eta1 eta2 eta3'
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
    prop_values = '0.5 300e-9 5e-10 1.5 16.29e-6 0' #J/m^2 m - ?
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
    #outputs = exodus_out
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
    outputs = exodus_out
    output_properties = mu
  [../]
  [./L_cu_eps] #TODO: lägg in args
    type = ParsedMaterial
    material_property_names = 'mu kappa D_cu D_eps A_cu A_eps c_cu_eps c_eps_cu length_scale energy_scale time_scale'
    f_name = L_cu_eps
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_cu/A_cu+D_eps/A_eps)/(3*kappa*(c_cu_eps-c_eps_cu)^2)' #nm^3/eVs
  [../]
  [./L_cu_eta]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_cu D_eta A_cu A_eta c_cu_eta c_eta_cu length_scale energy_scale time_scale'
    f_name = L_cu_eta
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_cu/A_cu+D_eta/A_eta)/(3*kappa*(c_cu_eta-c_eta_cu)^2)' #nm^3/eVs
  [../]
  [./L_cu_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_cu A_sn A_cu c_sn_cu c_cu_sn length_scale energy_scale time_scale'
    f_name = L_cu_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_cu/A_cu)/(3*kappa*(c_sn_cu-c_cu_sn)^2)'
  [../]
  [./L_eps_eta]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_eps D_eta A_eps A_eta c_eps_eta c_eta_eps length_scale energy_scale time_scale'
    f_name = L_eps_eta
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_eps/A_eps+D_eta/A_eta)/(3*kappa*(c_eps_eta-c_eta_eps)^2)'
    #function = '0'
  [../]
  [./L_eps_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_eps D_sn A_eps A_sn c_eps_sn c_sn_eps length_scale energy_scale time_scale'
    f_name = L_eps_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_eps/A_eps+D_sn/A_sn)/(3*kappa*(c_eps_sn-c_sn_eps)^2)'
  [../]
  [./L_eta_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_eta A_sn A_eta c_sn_eta c_eta_sn length_scale energy_scale time_scale'
    f_name = L_eta_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_eta/A_eta)/(3*kappa*(c_sn_eta-c_eta_sn)^2)'
  [../]

  #[./L_imc_imc]
  #  type = ParsedMaterial
  #  material_property_names = 'L_cu_imc L_imc_sn'
  #  f_name = L_imc_imc
  #  function = 'L_cu_imc'
  #[../]
  #Free energy
  [./fch_cu] #Chemical energy Cu phase
      type = DerivativeParsedMaterial
      f_name = fch0
      args = 'c0'
      material_property_names = 'A_cu B_cu C_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c0-chat_cu)^2+B_cu*(c0-chat_cu)+C_cu)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_imc1] #Chemical energy Cu3Sn phase
      type = DerivativeParsedMaterial
      f_name = fch1
      args = 'c1'
      material_property_names = 'A_eps B_eps C_eps chat_eps length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eps*(c1-chat_eps)^2+B_eps*(c1-chat_eps)+C_eps)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_imc2] #Chemical energy Cu6Sn5 phase grain 2
      type = DerivativeParsedMaterial
      f_name = fch2
      args = 'c2'
      material_property_names = 'A_eta B_eta C_eta chat_eta length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eta*(c2-chat_eta)^2+B_eta*(c2-chat_eta)+C_eta)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_sn] #Chemical energy Sn phase
      type = DerivativeParsedMaterial
      f_name = fch3
      args = 'c3'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c3-chat_sn)^2+B_sn*(c3-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
  [../]
  #SwitchingFunction
  [./h_cu]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta0
      outputs = exodus_out
      output_properties = h0
  [../]

  [./h_imc1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta1
      outputs = exodus_out
      output_properties = h1
  [../]
  [./h_imc2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta2
      outputs = exodus_out
      output_properties = h2
  [../]

  [./h_sn]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta3
      outputs = exodus_out
      output_properties = h3
  [../]

  #Double well, not used MAYBE USE TO KEEP THE ORDER PARAMETERS IN [0:1]
  [./g0]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta0
    well_only = True
    function_name = g0
  [../]
  #Double well, not used
  [./g1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta1
    well_only = True
    function_name = g1
  [../]

  [./g2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta2
    well_only = True
    function_name = g2
  [../]
  #Double well, not used
  [./g3]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta=eta3
    well_only = True
    function_name = g3
  [../]

  [./Mgb]
    type=ParsedMaterial
    material_property_names = 'D_gb delta delta_real h0(eta0,eta1,eta2,eta3) h1(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3) h3(eta0,eta1,eta2,eta3) A_cu A_eps A_eta A_sn length_scale energy_scale time_scale'
    f_name = Mgb
    function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h1*A_eps+h2*A_eta+h3*A_sn)*delta)'
    #function = '4e-5'
  [../]
  [./CHMobility]
      type = DerivativeParsedMaterial
      f_name = M
      args = 'eta0 eta1 eta2 eta3'
      material_property_names = 'h0(eta0,eta1,eta2,eta3) h1(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3) h3(eta0,eta1,eta2,eta3) D_cu D_eps D_eta D_sn A_cu A_eps A_eta A_sn Mgb length_scale energy_scale time_scale'
      #function = 's:=eta_cu^2+eta_imc1^2+eta_imc2^2+eta_sn^2;p:=eta_imc1^2*eta_imc2^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc1*h_imc2>1./16.,0,Mgb)' #nm^5/eVs
      function = '(length_scale^5/(energy_scale*time_scale))*(h0*D_cu/A_cu+h1*D_eps/A_eps+h2*D_eta/A_eta+h3*D_sn/A_sn)' #'+h_imc1*h_imc2*Mgb' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
      derivative_order = 2
      #outputs = exodus_out
  [../]

  [./ACMobility]
      type = DerivativeParsedMaterial
      f_name = L
      args = 'eta0 eta1 eta2 eta3'
      #material_property_names = 'L_cu_eps L_cu_eta L_cu_sn L_eps_eta L_eps_sn L_eta_sn'
      material_property_names = 'L_cu_eps L_cu_eta L_cu_sn L_eps_eta L_eps_sn L_eta_sn'
      # Added epsilon to prevent division by 0 (Larry Aagesen)
      #function ='pf:=1e5;eps:=0.01;(L_cu_eps*(pf*eta_cu^2+eps)*(pf*eta_imc1^2+eps)+L_cu_eta*(pf*eta_cu^2+eps)*(pf*eta_imc2^2+eps)+L_eps_sn*(pf*eta_imc1^2+eps)*(pf*eta_sn^2+eps)+L_eps_eta*(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps)+L_eta_sn*(pf*eta_imc2^2+eps)*(pf*eta_sn^2+eps)+L_cu_sn*(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps))/((pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))'
      function ='pf:=1e5;eps:=1e-5;(L_cu_eps*(pf*eta0^2+eps)*(pf*eta1^2+eps)+L_cu_eta*(pf*eta0^2+eps)*(pf*eta2^2+eps)+L_cu_sn*(pf*eta0^2+eps)*(pf*eta3^2+eps)+L_eps_eta*(pf*eta1^2+eps)*(pf*eta2^2+eps)+L_eps_sn*(pf*eta1^2+eps)*(pf*eta3^2+eps)+L_eta_sn*(pf*eta2^2+eps)*(pf*eta3^2+eps))/((pf*eta0^2+eps)*(pf*eta1^2+eps)+(pf*eta0^2+eps)*(pf*eta2^2+eps)+(pf*eta0^2+eps)*(pf*eta3^2+eps)+(pf*eta1^2+eps)*(pf*eta2^2+eps)+(pf*eta1^2+eps)*(pf*eta3^2+eps)+(pf*eta2^2+eps)*(pf*eta3^2+eps))'
      #function ='L_imc_sn'

      # Conditional function (Daniel Schwen)
      #function ='numer:=L_cu_eps*eta0^2*eta1^2+L_eps_eta*eta1^2*eta2^2+L_eta_sn*eta2^2*eta3^2;denom:=eta0^2*eta1^2+eta1^2*eta2^2+eta2^2*eta3^2;if(denom!=0,numer/denom,0.)'
      #function ='numer:=L_cu_eps*eta0^2*eta1^2+L_eps_eta*eta1^2*eta2^2+L_eta_sn*eta2^2*eta3^2;denom:=eta0^2*eta1^2+eta1^2*eta2^2+eta2^2*eta3^2;if(denom>0.0001,numer/denom,0.)'

      derivative_order = 2
      #outputs = exodus_out
  [../]

  [./time]
    type = TimeStepMaterial
    prop_time = time
    prop_dt = dt
  [../]
  [./noise_constants]
    type = GenericConstantMaterial
    prop_names = 'T kb lambda dim' #temperature Boltzmann gridsize dimensionality
    prop_values = '493 8.6173303e-5 50 3'
  [../]
  [./nuc_eps]
    type =  DerivativeParsedMaterial
    f_name = nuc_eps
    #args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'time dt T kb lambda dim L h0 h2' #'h0(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3)'
    #function = '20*h0*h2'
    function = 'if(time>28800&time<28805&h0*h2>0.09,sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007) (without h0*h2)
    #function = '100'
    outputs = exodus_out
  [../]
  [./nuc_eta]
    type =  DerivativeParsedMaterial
    f_name = nuc_eta
    #args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'time dt T kb lambda dim L h0 h3' #'h0(eta0,eta1,eta2,eta3) h2(eta0,eta1,eta2,eta3)'
    #function = '20*h0*h2'
    #function = 'sqrt(2*kb*T*L/(lambda^dim*dt))*h0*h3' #expression from Shen (2007) (without h0*h2)
    function = 'if(h0*h3>0.09&time<1,sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007) (without h0*h2)
    #function = '100'
    outputs = exodus_out
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
    displacements = 'disp_x disp_y disp_z'
    eigenstrain_names = eT_eta
  [../]
  [./stress_eta]
    type = ComputeLinearElasticStress
    base_name = eta2
  [../]
  [./eigenstrain_eta]
    type = ComputeEigenstrain
    base_name = eta2
    eigen_base = '1 1 1'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eta
    prefactor = pre
  [../]
  [./pre]
    type = ParsedMaterial
    args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'time h2(eta0,eta1,eta2,eta3)'
    f_name = pre
    #function = 'if(time<1,-0.003*time,-0.003)'
    #function = '-0.003'
    function = '0.02'
    #function = '0'
  [../]
  [./fel_eta]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta2
    f_name = fel2
    outputs = exodus_out
    output_properties = fel2
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
    displacements = 'disp_x disp_y disp_z'
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
    outputs = exodus_out
    output_properties = fel1
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
  [./fel_cu]
    type = ElasticEnergyMaterial
    args = ' '
    base_name = eta0
    f_name = fel0
    outputs = exodus_out
    output_properties = fel0
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
    displacements = 'disp_x disp_y disp_z'
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
    outputs = exodus_out
    output_properties = fel3
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
    displacements = 'disp_x disp_y disp_z'
  [../]

  # Total elastic energy
  [./Total_elastic_energy]
    type = DerivativeMultiPhaseMaterial
    f_name = fel_tot
    W = 0
    fi_names = 'fel0 fel1 fel2 fel3'
    hi_names = 'h0 h1 h2 h3'
    etas = 'eta0 eta1 eta2 eta3'
    outputs = exodus_out
    output_properties = fel_tot
    g = g0
  [../]


  #sum chemical and elastic energies
  [./F_cu]
    type = DerivativeSumMaterial
    f_name = F0
    args = 'c0'
    sum_materials = 'fch0 fel0'
    #sum_materials = 'fch0'
  [../]
  [./F_eps]
    type = DerivativeSumMaterial
    f_name = F1
    args = 'c1'
    sum_materials = 'fch1 fel1'
    #sum_materials = 'fch1'
  [../]
  [./F_eta]
    type = DerivativeSumMaterial
    f_name = F2
    args = 'c2'
    sum_materials = 'fch2 fel2'
    #sum_materials = 'fch2'
  [../]
  [./F_sn]
    type = DerivativeSumMaterial
    f_name = F3
    args = 'c3'
    sum_materials = 'fch3 fel3'
    #sum_materials = 'fch3'
  [../]

[]

[Kernels]
    [./TensorMechanics]
      displacements = 'disp_x disp_y disp_z'
      base_name = global
    [../]
    #Kernels for split Cahn-Hilliard equation
    # Cahn-Hilliard Equation
    [./CHBulk] # Gives the residual for the concentration, dF/dc-mu
        type = KKSSplitCHCRes
        variable = c
        ca       = c2
        cb       = c3
        fa_name  = F2 #only fa is used
        fb_name  = F3
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
        args = 'eta0 eta1 eta2 eta3'
    [../]

    #KKS conditions
    # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
    [./chempot_cu_imc]
      type = KKSPhaseChemicalPotential
      variable = c0
      cb       = c1
      fa_name  = F0
      fb_name  = F1
    [../]
    [./chempot_imc_imc]
      type = KKSPhaseChemicalPotential
      variable = c1
      cb       = c2
      fa_name  = F1
      fb_name  = F2
    [../]
    [./chempot_sn_cu]
      type = KKSPhaseChemicalPotential
      variable = c2
      cb       = c3
      fa_name  = F2
      fb_name  = F3
    [../]
    [./phaseconcentration] # enforce c = sum h_i*c_i
      type = KKSMultiPhaseConcentration
      variable = c3
      cj = 'c0 c1 c2 c3'
      hj_names = 'h0 h1 h2 h3'
      etas = 'eta0 eta1 eta2 eta3'
      c = c
    [../]

    #Kernels for Allen-Cahn equations
    [./KKSMultiACKernel]
      op_num = 4
      op_name_base = 'eta'
      ci_name_base = 'c'
      f_name_base = F
      #wi = 0.0624
      wi = 4.
      g_name_base = g
    [../]

    #Nucleation Kernel
    [./nucleation_eps]
      type = LangevinNoise # TODO: This draws random number from a uniform distribution, it should be from a standard gaussian instead(?)
      variable = eta1
      amplitude = 1
      seed = 123456789
      multiplier = nuc_eps
    [../]

    #Nucleation Kernel
    [./nucleation_eta]
      type = LangevinNoise
      variable = eta2
      amplitude = 1
      seed = 987654321
      multiplier = nuc_eta
    [../]
[]

[AuxVariables]
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
    [./s]
      order = FIRST
      family = LAGRANGE
    [../]
    [./von_mises]
      #Dependent variable used to visualize the Von Mises stress
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./sigma11]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./sigma22]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./sigma33]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./e11]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./e12]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./e22]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./e33]
      order = CONSTANT
      family = MONOMIAL
    [../]
[]

[AuxKernels]
  [./phase_id]
    type = ParsedAux
    variable = phase_id
    args = 'h0 h1 h2 h3'
    function = 'if(h0>0.5,1,0)+if(h1>0.5,2,0)+if(h2>0.5,3,0)+if(h3>0.5,4,0)'
    execute_on = timestep_end
  [../]

    [./von_mises_kernel]
      #Calculates the von mises stress and assigns it to von_mises
      type = RankTwoScalarAux
      variable = von_mises
      rank_two_tensor =global_stress
      execute_on = timestep_end
      scalar_type = VonMisesStress #TODO: Check units
    [../]
    [./matl_sigma11]
      type = RankTwoAux
      rank_two_tensor = global_stress
      index_i = 0
      index_j = 0
      variable = sigma11
    [../]
    [./matl_sigma22]
      type = RankTwoAux
      rank_two_tensor = global_stress
      index_i = 1
      index_j = 1
      variable = sigma22
    [../]
    [./matl_sigma33]
      type = RankTwoAux
      rank_two_tensor = global_stress
      index_i = 2
      index_j = 2
      variable = sigma33
    [../]
    [./matl_e11]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 0
      variable = e11
    [../]
    [./matl_e12]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 0
      index_j = 1
      variable = e12
    [../]
    [./matl_e22]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 1
      index_j = 1
      variable = e22
    [../]
    [./matl_e33]
      type = RankTwoAux
      rank_two_tensor = total_strain
      index_i = 2
      index_j = 2
      variable = e33
    [../]
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
        w = 4.
        execute_on = 'initial timestep_end'
    [../]
    [./f_int]
        type = ParsedAux
        variable = f_int
        args = 'eta0 eta1 eta2 eta3'
        constant_names = 'sigma delta gamma length_scale energy_scale'
        constant_expressions = '0.5 60e-9 1.5 1e9 6.24150943e18'
        function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta0^4-0.5*eta0^2+0.25*eta1^4-0.5*eta1^2+0.25*eta2^4-0.5*eta2^2+0.25*eta3^4-0.5*eta3^2+gamma*(eta0^2*(eta1^2+eta2^2+eta3^2)+eta1^2*(eta2^2+eta3^2)+eta2^2*eta3^2)+0.25)'
        execute_on = 'initial timestep_end'
    [../]
    [./s]
      type = ParsedAux
      variable = s
      args = 'eta0 eta1 eta2 eta3'
      #function = 'eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2'
      function = 'eta0+eta1+eta2+eta3'
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
      mat_prop = h0
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc1_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h1
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc2_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h2
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./sn_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h3
      execute_on = 'Initial TIMESTEP_END'
    [../]
    #Monitoring the progress
    #[./time]
    #  type = RunTime
    #  time_type = active
    #[../]
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
  line_search = none
  #petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  #petsc_options_value = 'asm       ilu            nonzero'
  #petsc_options_value = 'bjacobi       ilu            nonzero'

  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  l_max_its = 100
  nl_max_its = 20
  l_tol = 1.0e-4
  #l_abs_step_tol = 1e-8
  nl_rel_tol = 1.0e-7 #1.0e-10
  nl_abs_tol = 1.0e-8#1.0e-11

  #num_steps = 2000
  end_time = 180000 #50 hours
  #very simple adaptive time stepper
  scheme = bdf2
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
  #    # Turn on time stepping
  #    type =
  #    percent_change = 0.9
  #    dt = 0.002
  #
  #
  #[../]


[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    solve_type = PJFNK
    full = true
    #petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type  -sub_pc_factor_shift_type'
    #petsc_options_value = '     121              1e-10     1e-8     1e-5          lu   nonzero'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -sub_pc_type  -sub_pc_factor_shift_type'
    petsc_options_value = '     121              1e-10     1e-8     1e-5          lu   ilu    nonzero'
  [../]
  [./mydebug]
    type = FDP
    full = true
  [../]
[]

[Outputs]
  file_base = 3d-10x200x1-50nm-bothIMC-0.02-bdf2-w4
  print_linear_residuals = true
  [./exodus_out]
    type = Exodus
    interval = 1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  interval = 1 #5
[]
