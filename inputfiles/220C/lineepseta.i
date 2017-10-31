[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 1000
    ny = 1
    xmin = 0
    xmax = 100000 #[nm]
    ymin = 0
    ymax = 100
    elem_type = QUAD4
[]

[BCs]
    [./neumann]
        type = NeumannBC
        boundary = 'left right'
        variable = 'c w c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc1 eta_imc2 eta_sn'
        value = 0
    [../]
    [./Periodic]
      [./y]
        auto_direction = y
        variable = 'c w c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc1 eta_imc2 eta_sn'
      [../]
    [../]
[]

[Variables]
    # concentration Sn
    [./c]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]
    # chemical potential
    [./w]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]

    # phase concentration  Sn in Cu
    [./c_cu]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.02
        #initial_condition = 0.10569
        #scaling = 1e3
    [../]

    # phase concentration  Sn in Cu3Sn
    [./c_imc1]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.2433
        #scaling = 1e3
    [../]
    # phase concentration  Sn in Cu6Sn5
    [./c_imc2]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.4351
        #scaling = 1e3
    [../]

    # phase concentration  Sn in Sn
    [./c_sn]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.9889
        #scaling = 1e3
    [../]

    # order parameter Cu
    [./eta_cu]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]
    # order parameter Cu6Sn5
    [./eta_imc1]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]
    [./eta_imc2]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]

    # order parameter Sn
    [./eta_sn]
        order = FIRST
        family = LAGRANGE
        #scaling = 1e3
    [../]

[]
[ICs]
    [./eta_cu] #Cu
      type = FunctionIC
      variable = eta_cu
      function = 'if(x<=10000,1,0)'
    [../]
    [./eta_imc1] #Cu3Sn
      type = FunctionIC
      variable = eta_imc1
      function = 'if(x>10000&x<=11000,1,0)'
    [../]
    [./eta_imc2] #Cu6Sn5
      type = FunctionIC
      variable = eta_imc2
      function = 'if(x>11000&x<=12000,1,0)'
    [../]
    [./eta_sn] #Sn
      type = FunctionIC
      variable = eta_sn
      function = 'if(x>12000,1,0)'
    [../]

    [./c] #Concentration of Sn
      type = FunctionIC
      variable = c
      function = 'if(x<=10000,0.02,0)+if(x>10000&x<=11000,0.2433,0)+if(x>11000&x<=12000,0.4351,0)+if(x>12000,0.9889,0)'
      #function = 'c_cu*eta_cu+c_imc1*eta_imc1+c_imc2*eta_imc2+c_sn*eta_sn'
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
    prop_values = '0.5 400e-9 5e-10 1.5 16.29e-6 0' #J/m^2 m - ?
  [../]
  #Constants The energy parameters are for 220 C
  [./energy_A]
    type = GenericConstantMaterial
    prop_names = 'A_cu A_eps A_eta A_sn'
    #prop_values = '1.0133e5/Vm 4e5/Vm 4.2059e6/Vm' #J/m^3
    prop_values = '1.7756e10 2.4555e11 2.4555e11 2.3033e10' #J/m^3
  [../]
  [./energy_B]
    type = GenericConstantMaterial
    prop_names = 'B_cu B_eps B_eta B_sn'
    #prop_values = '-2.1146e4/Vm -6.9892e3/Vm 7.168e3/Vm' #J/m^3
    prop_values = '-2.6351e9 -1.4014e9 2.3251e7 2.14216e8' #J/m^3
  [../]
  [./energy_C]
    type = GenericConstantMaterial
    prop_names = 'C_cu C_eps C_eta C_sn'
    #prop_values = '-1.2842e4/Vm -1.9185e4/Vm -1.5265e4/Vm' #J/m^3
    prop_values = '-1.1441e9 -1.7294e9 -1.7646e9 -1.646e9' #J/m^3
  [../]
  [./energy_c_ab]
    type = GenericConstantMaterial
    prop_names = 'c_cu_eps c_cu_eta c_cu_sn c_eps_cu c_eps_eta c_eps_sn c_eta_cu c_eta_eps c_eta_sn c_sn_cu c_sn_eps c_sn_eta'
    prop_values = '0.02 0.1957 0.6088 0.2383 0.2483 0.000 0.4299 0.4343 0.4359 0.9789 0.000 0.9889' #-
  [../]
  [./energy_chat]
    type = GenericConstantMaterial
    prop_names = 'chat_cu chat_eps chat_eta chat_sn'
    prop_values = '0.02 0.2433 0.4351 0.9889' #-
  [../]
  [./diffusion_constants]
    type = GenericConstantMaterial
    prop_names = 'D_cu D_eps D_eta D_sn'
    #prop_values = '2.877e-36 6.575e-19 2.452e-17' # m^2/s
    prop_values = '2.8311e-23 2.7650e-15 5.9702e-15 4.0866e-14' # m^2/s
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
  [./L_eps_eta]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_eps D_eta A_eps A_eta c_eps_eta c_eta_eps length_scale energy_scale time_scale'
    f_name = L_eps_eta
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_eps/A_eps+D_eta/A_eta)/(3*kappa*(c_eps_eta-c_eta_eps)^2)'
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
  [./L_cu_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_cu A_sn A_cu c_sn_cu c_cu_sn length_scale energy_scale time_scale'
    f_name = L_cu_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_cu/A_cu)/(3*kappa*(c_sn_cu-c_cu_sn)^2)'
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
      f_name = fch_cu
      args = 'c_cu'
      material_property_names = 'A_cu B_cu C_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c_cu-chat_cu)^2+B_cu*(c_cu-chat_cu)+C_cu)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_imc1] #Chemical energy Cu3Sn phase
      type = DerivativeParsedMaterial
      f_name = fch_imc1
      args = 'c_imc1'
      material_property_names = 'A_eps B_eps C_eps chat_eps length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eps*(c_imc1-chat_eps)^2+B_eps*(c_imc1-chat_eps)+C_eps)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_imc2] #Chemical energy Cu6Sn5 phase grain 2
      type = DerivativeParsedMaterial
      f_name = fch_imc2
      args = 'c_imc2'
      material_property_names = 'A_eta B_eta C_eta chat_eta length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eta*(c_imc2-chat_eta)^2+B_eta*(c_imc2-chat_eta)+C_eta)' #eV/nm^3
      derivative_order = 2
  [../]
  [./fch_sn] #Chemical energy Sn phase
      type = DerivativeParsedMaterial
      f_name = fch_sn
      args = 'c_sn'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c_sn-chat_sn)^2+B_sn*(c_sn-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
  [../]
    #SwitchingFunction
    [./h_cu]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h_cu
        all_etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        phase_etas = eta_cu
    [../]

    [./h_imc1]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h_imc1
        all_etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        phase_etas = eta_imc1
    [../]
    [./h_imc2]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h_imc2
        all_etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        phase_etas = eta_imc2
    [../]

    [./h_sn]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h_sn
        all_etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        phase_etas = eta_sn
    [../]

    #Double well, not used MAYBE USE TO KEEP THE ORDER PARAMETERS IN [0:1]
    [./g_cu]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_cu
      well_only = True
      function_name = g_cu
    [../]
    #Double well, not used
    [./g_imc1]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_imc1
      well_only = True
      function_name = g_imc1
    [../]
    [./g_imc2]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_imc2
      well_only = True
      function_name = g_imc2
    [../]
    #Double well, not used
    [./g_sn]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_sn
      well_only = True
      function_name = g_sn
    [../]
    [./Mgb]
      type=ParsedMaterial
      material_property_names = 'D_gb delta delta_real h_cu(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc1(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc2(eta_cu,eta_imc1,eta_imc2,eta_sn) h_sn(eta_cu,eta_imc1,eta_imc2,eta_sn) A_cu A_eps A_eta A_sn length_scale energy_scale time_scale'
      f_name = Mgb
      function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h_cu*A_cu+h_imc1*A_eps+h_imc2*A_eta+h_sn*A_sn)*delta)'
      #function = '4e-5'
    [../]
    [./CHMobility]
        type = DerivativeParsedMaterial
        f_name = M
        args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        material_property_names = 'h_cu(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc1(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc2(eta_cu,eta_imc1,eta_imc2,eta_sn) h_sn(eta_cu,eta_imc1,eta_imc2,eta_sn) D_cu D_eps D_eta D_sn A_cu A_eps A_eta A_sn Mgb length_scale energy_scale time_scale'
        #function = 's:=eta_cu^2+eta_imc1^2+eta_imc2^2+eta_sn^2;p:=eta_imc1^2*eta_imc2^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
        #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc1*h_imc2>1./16.,0,Mgb)' #nm^5/eVs
        function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_eps/A_eps+h_imc2*D_eta/A_eta+h_sn*D_sn/A_sn)' #'+h_imc1*h_imc2*Mgb' #nm^5/eVs
        #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
        derivative_order = 2
        outputs = exodus_out
    [../]

    [./ACMobility]
        type = DerivativeParsedMaterial
        f_name = L
        args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        #material_property_names = 'L_cu_eps L_cu_eta L_cu_sn L_eps_eta L_eps_sn L_eta_sn'
        material_property_names = 'L_cu_eps L_eps_eta L_eta_sn'
        # Added epsilon to prevent division by 0 (Larry Aagesen)
        #function ='pf:=1e5;eps:=0.01;(L_cu_eps*(pf*eta_cu^2+eps)*(pf*eta_imc1^2+eps)+L_cu_eta*(pf*eta_cu^2+eps)*(pf*eta_imc2^2+eps)+L_eps_sn*(pf*eta_imc1^2+eps)*(pf*eta_sn^2+eps)+L_eps_eta*(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps)+L_eta_sn*(pf*eta_imc2^2+eps)*(pf*eta_sn^2+eps)+L_cu_sn*(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps))/((pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))'
        function ='pf:=1e5;eps:=0.01;(L_cu_eps*(pf*eta_cu^2+eps)*(pf*eta_imc1^2+eps)+L_eps_eta*(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps)+L_eta_sn*(pf*eta_imc2^2+eps)*(pf*eta_sn^2+eps))/((pf*eta_cu^2+eps)*(pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))'
        #function ='L_imc_sn'

        # Conditional function (Daniel Schwen)
        #function ='numer:=L_cu_imc*eta_cu^2*eta_imc^2+L_imc_sn*eta_imc^2*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2;denom:=eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2;if(denom!=0,numer/denom,0.5*(L_cu_imc+L_imc_sn))'

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
        ca       = c_imc1
        cb       = c_sn
        fa_name  = fch_imc1 #only fa is used
        fb_name  = fch_sn
        w        = w
        h_name   = h_imc1
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
        args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
    [../]

    #KKS conditions
    # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
    [./chempot_cu_imc]
      type = KKSPhaseChemicalPotential
      variable = c_cu
      cb       = c_imc1
      fa_name  = fch_cu
      fb_name  = fch_imc1
    [../]
    [./chempot_imc_imc]
      type = KKSPhaseChemicalPotential
      variable = c_imc1
      cb       = c_imc2
      fa_name  = fch_imc1
      fb_name  = fch_imc2
    [../]
    [./chempot_sn_cu]
      type = KKSPhaseChemicalPotential
      variable = c_imc2
      cb       = c_sn
      fa_name  = fch_imc2
      fb_name  = fch_sn
    [../]
    [./phaseconcentration] # enforce c = sum h_i*c_i
      type = KKSMultiPhaseConcentration
      variable = c_sn
      cj = 'c_cu c_imc1 c_imc2 c_sn'
      hj_names = 'h_cu h_imc1 h_imc2 h_sn'
      etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
      c = c
    [../]

    #Kernels for Allen-Cahn equation for Cu
    [./detadt_cu]
      type = TimeDerivative
      variable = eta_cu
    [../]
    [./ACBulkF_cu] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used(?)
      type = KKSMultiACBulkF
      variable  = eta_cu
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      gi_name   = g_cu
      eta_i     = eta_cu
      wi        = 10
      mob_name = L
      args      = 'c_cu c_imc1 c_imc2 c_sn eta_imc1 eta_imc2 eta_sn'
    [../]
    [./ACBulkC_cu] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_cu
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      cj_names  = 'c_cu c_imc1 c_imc2 c_sn'
      eta_i     = eta_cu
      mob_name = L
      args      = 'eta_imc1 eta_imc2 eta_sn'
    [../]
    [./ACInterface_cu] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_cu
      kappa_name = kappa
      mob_name = L
      args      = 'eta_imc1 eta_imc2 eta_sn'
      variable_L = true
    [../]
    [./ACdfintdeta_cu] #L*m*(eta_i^3-eta_i+2*beta*eta_i*sum_j eta_j^2)
      type = ACGrGrMulti
      variable = eta_cu
      v = 'eta_imc1 eta_imc2 eta_sn'
      gamma_names = 'gamma gamma gamma'
      mob_name = L
      args = 'eta_imc1 eta_imc2 eta_sn'
    [../]

    #Kernels for Allen-Cahn equation for Cu6Sn5
    [./detadt_imc1]
      type = TimeDerivative
      variable = eta_imc1
    [../]
    [./ACBulkF_imc1] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_imc1
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      gi_name   = g_imc1
      eta_i     = eta_imc1
      wi        = 10
      mob_name = L
      args      = 'c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc2 eta_sn'
    [../]
    [./ACBulkC_imc1] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_imc1
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      cj_names  = 'c_cu c_imc1 c_imc2 c_sn'
      eta_i     = eta_imc1
      mob_name = L
      args      = 'eta_cu eta_imc2 eta_sn'
    [../]
    [./ACInterface_imc1] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_imc1
      kappa_name = kappa
      mob_name = L
      args      = 'eta_cu eta_imc2 eta_sn'
      variable_L = true
    [../]
    [./ACdfintdeta_imc1]
      type = ACGrGrMulti
      variable = eta_imc1
      v = 'eta_cu eta_imc2 eta_sn'
      gamma_names = 'gamma gamma gamma'
      mob_name = L
      args = 'eta_cu eta_imc2 eta_sn'
    [../]

    [./detadt_imc2]
      type = TimeDerivative
      variable = eta_imc2
    [../]
    [./ACBulkF_imc2] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_imc2
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      gi_name   = g_imc2
      eta_i     = eta_imc2
      wi        = 10
      mob_name = L
      args      = 'c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc1 eta_sn'
    [../]
    [./ACBulkC_imc2] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_imc2
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      cj_names  = 'c_cu c_imc1 c_imc2 c_sn'
      eta_i     = eta_imc2
      mob_name = L
      args      = 'eta_cu eta_imc1 eta_sn'
    [../]
    [./ACInterface_imc2] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_imc2
      kappa_name = kappa
      mob_name = L
      args      = 'eta_cu eta_imc1 eta_sn'
      variable_L = true
    [../]
    [./ACdfintdeta_imc2]
      type = ACGrGrMulti
      variable = eta_imc2
      v = 'eta_cu eta_imc1 eta_sn'
      gamma_names = 'gamma gamma gamma'
      mob_name = L
      args = 'eta_cu eta_imc1 eta_sn'
    [../]
    #Kernels for Allen-Cahn equation for Sn
    [./detadt_sn]
      type = TimeDerivative
      variable = eta_sn
    [../]
    [./ACBulkF_sn] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_sn
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      gi_name   = g_sn
      eta_i     = eta_sn
      wi        = 10
      mob_name = L
      args      = 'c_cu c_imc1 c_imc2 c_sn eta_imc1 eta_imc2 eta_cu'
    [../]
    [./ACBulkC_sn] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_sn
      Fj_names  = 'fch_cu fch_imc1 fch_imc2 fch_sn'
      hj_names  = 'h_cu h_imc1 h_imc2 h_sn'
      cj_names  = 'c_cu c_imc1 c_imc2 c_sn'
      eta_i     = eta_sn
      mob_name = L
      args      = 'eta_cu eta_imc1 eta_imc2'
    [../]
    [./ACInterface_sn] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_sn
      kappa_name = kappa
      mob_name = L
      args      = 'eta_cu eta_imc1 eta_imc2'
      variable_L = true
    [../]
    [./ACdfintdeta_sn]
      type = ACGrGrMulti
      variable = eta_sn
      v = 'eta_cu eta_imc1 eta_imc2'
      gamma_names = 'gamma gamma gamma'
      mob_name = L
      args= 'eta_cu eta_imc1 eta_imc2'
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
        hj_names = 'h_cu h_imc1 h_imc2 h_sn'
        Fj_names = 'fch_cu fch_imc1 fch_imc2 fch_sn'
        gj_names = 'g_cu g_imc1 g_imc2 g_sn'
        additional_free_energy = f_int
        interfacial_vars = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        kappa_names = 'kappa kappa kappa kappa'
        w = 10
        execute_on = 'initial timestep_end'
    [../]
    [./f_int]
        type = ParsedAux
        variable = f_int
        args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
        constant_names = 'sigma delta gamma length_scale energy_scale'
        constant_expressions = '0.5 0.667e-6 1.5 1e9 6.24150943e18'
        function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta_cu^4-0.5*eta_cu^2+0.25*eta_imc1^4-0.5*eta_imc1^2+0.25*eta_imc2^4-0.5*eta_imc2^2+0.25*eta_sn^4-0.5*eta_sn^2+gamma*(eta_cu^2*(eta_imc1^2+eta_imc2^2+eta_sn^2)+eta_imc1^2*(eta_imc2^2+eta_sn^2))+0.25)'
        execute_on = 'initial timestep_end'
    [../]
    [./s]
      type = ParsedAux
      variable = s
      args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
      #function = 'eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2'
      function = 'eta_cu+eta_imc1+eta_imc2+eta_sn'
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
      mat_prop = h_cu
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc1_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h_imc1
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./imc2_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h_imc2
      execute_on = 'Initial TIMESTEP_END'
    [../]
    [./sn_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h_sn
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
  end_time = 172800 #48 hours
  #very simple adaptive time stepper
  [./TimeStepper]
      # Turn on time stepping
      type = IterationAdaptiveDT
      dt = 1e-4
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
  file_base = line220C
  [./exodus_out]
    type = Exodus
    interval = 1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  interval = 1 #5
[]
