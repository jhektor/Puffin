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
        variable = 'c w c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc1 eta_imc2 eta_sn'
        value = 0
    [../]
    [./Periodic]
      [./x]
        auto_direction = 'x'
        variable = 'c w c_cu c_imc1 c_imc2 c_sn eta_cu eta_imc1 eta_imc2 eta_sn'
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
    [./c_cu]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.002
        #initial_condition = 0.10569
        #scaling = 1e2
    [../]

    # phase concentration  Sn in Cu6Sn5
    [./c_imc1]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.417
        #scaling = 1e2
    [../]
    # phase concentration  Sn in Cu6Sn5
    [./c_imc2]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.417
        #scaling = 1e2
    [../]

    # phase concentration  Sn in Sn
    [./c_sn]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.999
        #scaling = 1e2
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
      type = UnitySubVarIC
      variable = eta_cu
      etas = 'eta_imc1 eta_imc2 eta_sn'
    [../]
    [./eta_imc1] #Cu6Sn5
      type = SmoothSuperellipsoidIC
      variable = eta_imc1
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
      variable = eta_imc2
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
      variable = eta_sn
      etas = 'eta_imc1 eta_imc2'
      use_threshold = true
      y_threshold = 800
    [../]

    [./c] #Concentration of Sn
      type = VarDepIC
      variable = c
      etas = 'eta_cu eta_imc1 eta_imc2 eta_sn'
      cis = 'c_cu c_imc1 c_imc2 c_sn'
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
    prop_names = 'c_cu_imc c_imc_cu c_imc_sn c_sn_imc c_cu_sn c_sn_cu'
    prop_values = '0.1057 0.3821 0.4529 0.9994 0.3112 0.9976' #-
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
  [./L_cu_imc]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_cu D_imc A_cu A_imc c_cu_imc c_imc_cu length_scale energy_scale time_scale'
    f_name = L_cu_imc
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_cu/A_cu+D_imc/A_imc)/(3*kappa*(c_cu_imc-c_imc_cu)^2)' #nm^3/eVs
  [../]
  [./L_imc_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_imc A_sn A_imc c_sn_imc c_imc_sn length_scale energy_scale time_scale'
    f_name = L_imc_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_imc/A_imc)/(3*kappa*(c_sn_imc-c_imc_sn)^2)'
  [../]
  [./L_cu_sn]
    type = ParsedMaterial
    material_property_names = 'mu kappa D_sn D_cu A_sn A_cu c_sn_cu c_cu_sn length_scale energy_scale time_scale'
    f_name = L_cu_sn
    function = '(length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_cu/A_cu)/(3*kappa*(c_sn_cu-c_cu_sn)^2)'
  [../]
  [./L_imc_imc]
    type = ParsedMaterial
    material_property_names = 'L_cu_imc L_imc_sn'
    f_name = L_imc_imc
    function = 'L_cu_imc'
  [../]
  #Free energy
  [./fch_cu] #Chemical energy Cu phase
      type = DerivativeParsedMaterial
      f_name = fch_cu
      args = 'c_cu'
      material_property_names = 'A_cu B_cu Cf_cu chat_cu length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c_cu-chat_cu)^2+B_cu*(c_cu-chat_cu)+Cf_cu)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch_imc1] #Chemical energy Cu6Sn5 phase grain 1
      type = DerivativeParsedMaterial
      f_name = fch_imc1
      args = 'c_imc1'
      material_property_names = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c_imc1-chat_imc)^2+B_imc*(c_imc1-chat_imc)+Cf_imc)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch_imc2] #Chemical energy Cu6Sn5 phase grain 2
      type = DerivativeParsedMaterial
      f_name = fch_imc2
      args = 'c_imc2'
      material_property_names = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c_imc2-chat_imc)^2+B_imc*(c_imc2-chat_imc)+Cf_imc)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
  [../]
  [./fch_sn] #Chemical energy Sn phase
      type = DerivativeParsedMaterial
      f_name = fch_sn
      args = 'c_sn'
      material_property_names = 'A_sn B_sn Cf_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c_sn-chat_sn)^2+B_sn*(c_sn-chat_sn)+Cf_sn)' #eV/nm^3
      derivative_order = 2
      outputs = exodus_out
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
    #args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
    material_property_names = 'D_gb delta delta_real h_cu(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc1(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc2(eta_cu,eta_imc1,eta_imc2,eta_sn) h_sn(eta_cu,eta_imc1,eta_imc2,eta_sn) A_cu A_imc A_sn length_scale energy_scale time_scale'
    f_name = Mgb
    function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h_cu*A_cu+h_imc1*A_imc+h_imc2*A_imc+h_sn*A_sn)*delta)'
    #function = '4e-5'
  [../]
  [./CHMobility]
    type = DerivativeParsedMaterial
    f_name = M
    args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
    material_property_names = 'h_cu(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc1(eta_cu,eta_imc1,eta_imc2,eta_sn) h_imc2(eta_cu,eta_imc1,eta_imc2,eta_sn) h_sn(eta_cu,eta_imc1,eta_imc2,eta_sn) D_cu D_imc D_sn A_cu A_imc A_sn length_scale energy_scale time_scale Mgb'
    #function = 's:=eta_cu^2+eta_imc1^2+eta_imc2^2+eta_sn^2;p:=eta_imc1^2*eta_imc2^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc1*h_imc2>1./16.,0,Mgb)' #nm^5/eVs
    function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc1*D_imc/A_imc+h_imc2*D_imc/A_imc+h_sn*D_sn/A_sn)' #'+h_imc1*h_imc2*(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h_cu*A_cu+h_imc1*A_imc+h_imc2*A_imc+h_sn*A_sn)*delta)' #nm^5/eVs
    #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
    derivative_order = 2
    outputs = exodus_out
  [../]

  [./ACMobility]
    type = DerivativeParsedMaterial
    f_name = L
    args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
    material_property_names = 'L_cu_imc L_imc_sn L_cu_sn L_imc_imc' # h_cu(eta_cu,eta_imc,eta_sn) h_imc(eta_cu,eta_imc,eta_sn) h_sn(eta_cu,eta_imc,eta_sn)'

    # Added epsilon to prevent division by 0 (Larry Aagesen)
    #function ='pf:=1e5;eps:=0.01;(L_cu_imc*(pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+L_imc_sn*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+L_cu_sn*(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+L_imc_imc*(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))/((pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))'

    #function ='pf:=1e5;eps:=1e-5;(L_cu_imc*(pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+L_imc_sn*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+L_cu_sn*(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+L_imc_imc*(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))/((pf*eta_cu^2+eps)*((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))+((pf*eta_imc1^2+eps)+(pf*eta_imc2^2+eps))*(pf*eta_sn^2+eps)+(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_imc1^2+eps)*(pf*eta_imc2^2+eps))'
    function ='L_imc_sn'

    # Conditional function (Daniel Schwen)
    #function ='numer:=L_cu_imc*eta_cu^2*(eta_imc1^2+eta_imc2^2)+L_imc_sn*(eta_imc1^2+eta_imc2^2)*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2+L_imc_imc*eta_imc1^2*eta_imc2^2;denom:=eta_cu^2*(eta_imc1^2+eta_imc2^2)+(eta_imc1^2+eta_imc2^2)*eta_sn^2+eta_cu^2*eta_sn^2+eta_imc1^2*eta_imc^2;if(denom!=0,numer/denom,0.5*(L_cu_imc+L_imc_sn))'
    #function='numer:=L_cu_imc*eta_cu^2*(eta_imc1^2+eta_imc2^2)+L_imc_sn*(eta_imc1^2+eta_imc2^2)*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2+L_imc_imc*eta_imc1^2*eta_imc2^2;denom:=eta_cu^2*(eta_imc1^2+eta_imc2^2)+(eta_imc1^2+eta_imc2^2)*eta_sn^2+eta_cu^2*eta_sn^2+eta_imc1^2*eta_imc^2;if(denom!=0,numer/denom,0.)'
    #function='numer:=L_cu_imc*eta_cu^2*(eta_imc1^2+eta_imc2^2)+L_imc_sn*(eta_imc1^2+eta_imc2^2)*eta_sn^2+L_imc_imc*eta_imc1^2*eta_imc2^2;denom:=eta_cu^2*(eta_imc1^2+eta_imc2^2)+(eta_imc1^2+eta_imc2^2)*eta_sn^2+eta_imc1^2*eta_imc2^2;if(denom!=0,numer/denom,0.5*(L_cu_imc+L_imc_sn))' #0.0256 is 0.2^2*0.8^2
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
      wi        = 10.
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
      wi        = 10.
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
      wi        = 10.
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
      wi        = 10.
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
      w = 10.
      execute_on = 'initial timestep_end'
  [../]
  [./f_int]
      type = ParsedAux
      variable = f_int
      args = 'eta_cu eta_imc1 eta_imc2 eta_sn'
      constant_names = 'sigma delta gamma length_scale energy_scale'
      constant_expressions = '0.5 0.4e-6 1.5 1e9 6.24150943e18'
      function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta_cu^4-0.5*eta_cu^2+0.25*eta_imc1^4-0.5*eta_imc1^2+0.25*eta_imc2^4-0.5*eta_imc2^2+0.25*eta_sn^4-0.5*eta_sn^2+gamma*(eta_cu^2*(eta_imc1^2+eta_imc2^2+eta_sn^2)+eta_imc1^2*(eta_imc2^2+eta_sn^2)+eta_imc2^2*eta_sn^2)+0.25)'
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
  file_base = noaction-comp
  [./exodus_out]
    type = Exodus
    interval = 1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  interval = 1 #5
[]
