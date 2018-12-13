[Mesh]
  type = GeneratedMesh
  dim = 2
  elem_type = QUAD4 #HEX8
  nx = 200#100#80#80#12
  ny = 100#100#114#114#114#16
  # nz = 1
  xmin = -500#-3000
  xmax = 500#3000
  ymin = 0#-1000
  ymax = 500#7000
  # zmin = 0
  # zmax = 10
  displacements = 'disp_x disp_y' # disp_z'
[]
[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  displacements = 'disp_x disp_y' # disp_z'
[]

[Variables]
  [./disp_x]
    initial_condition = 0
    scaling = 1e-2
  [../]
  [./disp_y]
    initial_condition = 0
    scaling = 1e-2
  [../]
  # [./disp_z]
  #   initial_condition = 0
  #   scaling = 1e-2
  # [../]

  [./c]
    # scaling = 1e1
  [../]
  [./w]
    initial_condition = 0
  [../]

  #cu grain
  [./eta0]
  [../]
  [./c0]
    initial_condition = 0.002
    # scaling = 1e1
  [../]
  #sn grains
  #[100]
  [./eta1]
  [../]
  [./c1]
    initial_condition = 0.999
    # scaling = 1e1
  [../]
  #[110]
  [./eta2]
  [../]
  [./c2]
    initial_condition = 0.999
    # scaling = 1e1
  [../]
  #imc
  [./eta3]
  [../]
  [./c3]
    initial_condition = 0.417
    # scaling = 1e1
  [../]
[]
[ICs]
  [./eta0] # copper
    type = FunctionIC
    variable = eta0
    function = 'r:=sqrt(x^2+y^2); if(r<200,1,0)'
  [../]
  [./eta1] #sn1
    type = FunctionIC
    variable = eta1
    function = 'r:=sqrt(x^2+y^2); if(r>250 & x<=0,1,0)'
  [../]
  [./eta2] #sn2
    type = FunctionIC
    variable = eta2
    function = 'r:=sqrt(x^2+y^2); if(r>250 & x>0,1,0)'
  [../]
  [./eta3]
    type = UnitySubVarIC
    variable = eta3
    etas = 'eta0 eta1 eta2'
  [../]

  [./c]
    type = VarDepIC
    variable = c
    etas = 'eta0 eta1 eta2 eta3'
    cis = 'c0 c1 c2 c3'
    # cis = 'ccu csn csn csn cimc'
  [../]
[]

[BCs]
  # [./Periodic]
  #   [./xy]
  #     auto_direction = 'x z'
  #     variable = 'eta0 eta1 eta2 eta3 c w c0 c1 c2 c3 disp_x disp_y disp_z'
  #   [../]
  # [../]
  [./symy]
    type = PresetBC
    variable = disp_y
    value = 0
    # boundary = 'bottom'
    boundary = 'bottom top'
    # use_displaced_mesh = true
  [../]
  [./symx]
    type = PresetBC
    variable = disp_x
    value = 0
    boundary = 'left right'
    # boundary = 'bottom'
    # use_displaced_mesh = true
  [../]
  # [./symz]
  #   type = PresetBC
  #   variable = disp_z
  #   value = 0
  #   # boundary = 'bottom'
  #   boundary = 'front back'
  #   # use_displaced_mesh = true
  # [../]
[]

[UserObjects]
  #Crystal plasticity for the tin grains
  #State variables
  [./state_var1]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate1
    scale_factor = 1.0
  [../]
  [./state_var2]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate2
    scale_factor = 1.0
  [../]
  #State variable evolution rates
  [./state_var_evol_rate1]
    type = CrystalPlasticityStateVarRateComponentVoceP
    variable_size = 32
    groups = '0 2 4 6 10 20 24 32'
    h0_group_values = '0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484'
    tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    tauSat_group_values = '0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115'
    hardeningExponent_group_values = '2.0 2.0 2.0 2.0 2.0 2.0 2.0'
    coplanarHardening_group_values = '1.0 1.0 1.0 1.0 1.0 1.0 1.0' #q_aa = 1
    selfHardening_group_values = '1.4 1.4 1.4 1.4 1.4 1.4 1.4'
    crystal_lattice_type = BCT # default is BCT
    uo_slip_rate_name = slip_rate1
    uo_state_var_name = state_var1
  [../]
  [./state_var_evol_rate2]
    type = CrystalPlasticityStateVarRateComponentVoceP
    variable_size = 32
    groups = '0 2 4 6 10 20 24 32'
    h0_group_values = '0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484'
    tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    tauSat_group_values = '0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115'
    hardeningExponent_group_values = '2.0 2.0 2.0 2.0 2.0 2.0 2.0'
    coplanarHardening_group_values = '1.0 1.0 1.0 1.0 1.0 1.0 1.0' #q_aa = 1
    selfHardening_group_values = '1.4 1.4 1.4 1.4 1.4 1.4 1.4'
    crystal_lattice_type = BCT # default is BCT
    uo_slip_rate_name = slip_rate2
    uo_state_var_name = state_var2
  [../]
  #Slip rates
  [./slip_rate1]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct_miller.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    #flowprops = '1 32 0.001 2.6e-8' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var1
    base_name = eta1
  [../]
  [./slip_rate2]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct_miller.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    #flowprops = '1 32 0.001 2.6e-8' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var2
    base_name = eta2
  [../]
  # Slip resistance
  [./slip_resistance1]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var1
  [../]
  [./slip_resistance2]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var2
  [../]
[]

[Materials]   # PROPERTIES ARE FOR 150 DEGREE
  #Scalings
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale'
    prop_values = '1e9 6.24150943e18' #m->nm, J->eV
  [../]
  [./model_constants]
    type = GenericConstantMaterial # delta is diff reg. width
    prop_names = 'sigma delta delta_real gamma tgrad_corr_mult'
    prop_values = '0.5 50e-9 5e-10 1.5 0' #J/m^2 m - ?
    #prop_values = '0.5 60e-9 1.5 0' #J/m^2 m - ?
  [../]
  ## MECHANICS
  #Copper
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
  # Sn
  # [100] grain

  #Euler angles:
  # 0 0 0 should be c-axis out of screen
  # 0 90 90 should be a-axis out of screen

  [./elasticity_tensor1]
    type = ComputeElasticityTensorCPBaseName
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3' #MPa #From Darbandi 2014 table I
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    base_name = eta1
  [../]
  [./crysp1]
    type = FiniteStrainUObasedCPBaseName #hscale
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate1'
    uo_slip_resistances = 'slip_resistance1'
    uo_state_vars = 'state_var1'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate1'
    base_name = eta1
    maximum_substep_iteration = 10
    h_scale = h1
    h_tol = 0.001
  [../]
  [./strain1]
    type = ComputeFiniteStrain
    base_name = eta1
  [../]
  # [110] grain
  [./elasticity_tensor2]
    type = ComputeElasticityTensorCPBaseName
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3' #MPa #From Darbandi 2014 table I
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 90
    euler_angle_2 = 90
    euler_angle_3 = 0
    base_name = eta2
  [../]
  [./crysp2]
    type = FiniteStrainUObasedCPBaseName #hscale
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate2'
    uo_slip_resistances = 'slip_resistance2'
    uo_state_vars = 'state_var2'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate2'
    base_name = eta2
    maximum_substep_iteration = 10
    h_scale = h2
    h_tol = 0.001
  [../]
  [./strain2]
    type = ComputeFiniteStrain
    base_name = eta2
  [../]

  #Cu6Sn5
  [./elasticity_tensor_imc]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 701. #112.3 GPa in eV/nm^3
    poissons_ratio = 0.31
    base_name = eta3
  [../]
  [./strain_imc]
    type = ComputeFiniteStrain
    base_name = eta3
    eigenstrain_names = eigenstrain
  [../]
  [./stress_imc]
    type = ComputeFiniteStrainElasticStress
    base_name = eta3
  [../]
  [./eigenstrain]
    type = ComputeVariableEigenstrain
    args = 'eta0 eta1 eta2 eta3'
    base_name = eta3
    eigen_base = '1 1 1 0 0 0'
    eigenstrain_name = eigenstrain
    prefactor = prefactor
  [../]
  [./prefactor]
    type = DerivativeParsedMaterial
    args = 'imc_cu imc_sn'
    f_name = prefactor
    constant_names = 'esn ecu'
    constant_expressions = '-0.02 -0.02'
    # constant_expressions = '0.02 -0.02'
    function = 'esn*if(imc_sn>0.05,1,0)+ecu*if(imc_cu>0.05,1,0)'
    # function = '0'
    outputs = exodus
    output_properties = prefactor
  [../]
  #Global stress and strain
  [./global_stress] #homogeniserar bara Cauchy stress
    type = MultiPhaseStressMaterial
    phase_base = 'eta0 eta1 eta2 eta3'
    h          = 'h0 h1 h2 h3'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeFiniteStrain
    base_name = global
   # eigenstrain_names = eigenstrain
  [../]

  ## INTERPOLATION FUNCTIONS
  [./h0]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h0
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta0
      outputs = exodus
      output_properties = h0
      use_displaced_mesh = true
  [../]
  [./h1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta1
      outputs = exodus
      output_properties = h1
      use_displaced_mesh = true
  [../]
  [./h2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta2
      outputs = exodus
      output_properties = h2
      use_displaced_mesh = true
  [../]
  [./h3]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta3
      outputs = exodus
      output_properties = h3
      use_displaced_mesh = true
  [../]
  # Dummy double well
  [./g0]
    type = BarrierFunctionMaterial
    well_only = true
    function_name = g0
    g_order = SIMPLE
    eta = eta0
    use_displaced_mesh = true
  [../]
  [./g1]
    type = BarrierFunctionMaterial
    well_only = true
    function_name = g1
    g_order = SIMPLE
    eta = eta1
    use_displaced_mesh = true
  [../]
  [./g2]
    type = BarrierFunctionMaterial
    well_only = true
    function_name = g2
    g_order = SIMPLE
    eta = eta2
    use_displaced_mesh = true
  [../]
  [./g3]
    type = BarrierFunctionMaterial
    well_only = true
    function_name = g3
    g_order = SIMPLE
    eta = eta3
    use_displaced_mesh = true
  [../]

  ## FREE ENEERGYIES
  # Chemical free energies
  [./chemical_energy]
    type = GenericConstantMaterial
    prop_names = 'A_cu B_cu C_cu chat_cu
                  A_imc B_imc C_imc chat_imc
                  A_sn B_sn C_sn chat_sn
                  '
    prop_values = '5.8501e9 -1.378e9 -1.1749e9 0.144
                   2.4555e10 -5.3935e8 -1.5313e9 0.413
                   4.1436e10 2.9931e8 -1.3715e9 0.995'
  [../]
  [./fch_cu]  #Chemical energy Cu
    type = DerivativeParsedMaterial
    f_name = fch_cu
    args = c0
    material_property_names = 'A_cu B_cu C_cu chat_cu length_scale energy_scale'
    function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c0-chat_cu)^2+B_cu*(c0-chat_cu)+C_cu)'
    use_displaced_mesh = true
  [../]
  [./fch_sn100]  #Chemical energy Sn
    type = DerivativeParsedMaterial
    f_name = fch_sn1
    args = c1
    material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
    function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c1-chat_sn)^2+B_sn*(c1-chat_sn)+C_sn)'
    use_displaced_mesh = true
  [../]
  [./fch_sn110]  #Chemical energy Sn
    type = DerivativeParsedMaterial
    f_name = fch_sn2
    args = c2
    material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
    function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c2-chat_sn)^2+B_sn*(c2-chat_sn)+C_sn)'
    use_displaced_mesh = true
  [../]
  [./fch_imc]  #Chemical energy Sn
    type = DerivativeParsedMaterial
    f_name = fch_imc
    args = c3
    material_property_names = 'A_imc B_imc C_imc chat_imc length_scale energy_scale'
    function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c3-chat_imc)^2+B_imc*(c3-chat_imc)+C_imc)'
    use_displaced_mesh = true
  [../]
  # Elastic energies
  [./fe0]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe0
    base_name = eta0
    use_displaced_mesh = true
  [../]
  [./fe1]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe1
    base_name = eta1
    plasticity = true
    use_displaced_mesh = true
  [../]
  [./fe2]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe2
    base_name = eta2
    plasticity = true
    use_displaced_mesh = true
  [../]
  [./fe3]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe3
    base_name = eta3
    eigenstrain = true
    eigenstrain_name = eigenstrain
    use_displaced_mesh = true
  [../]
  # Plastic energies
  [./fp1]
    type = CPPlasticEnergyMaterial
    # Q = 1 in class
    variable_size = 32
    uo_state_var_name = state_var1
    f_name = fp1
    use_displaced_mesh = true
    args = ' '
    #s0 = 0.144
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    # outputs = exodus
    # output_properties = fp2
  [../]
  [./fp2]
    type = CPPlasticEnergyMaterial
    # Q = 1 in class
    variable_size = 32
    uo_state_var_name = state_var2
    f_name = fp2
    use_displaced_mesh = true
    args = ' '
    #s0 = 0.144
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    # outputs = exodus
    # output_properties = fp2
  [../]

  # Add energies together
  [./F0]
    type = DerivativeSumMaterial
    f_name = F0
    args = 'c0 eta0 eta1 eta2 eta3'
    sum_materials = 'fch_cu fe0'
    use_displaced_mesh = true
    # outputs = exodus
  [../]
  [./F1]
    type = DerivativeSumMaterial
    f_name = F1
    args = 'c1 eta0 eta1 eta2 eta3'
    sum_materials = 'fch_sn1 fe1 fp1'
    use_displaced_mesh = true
    # outputs = exodus
  [../]
  [./F2]
    type = DerivativeSumMaterial
    f_name = F2
    args = 'c2 eta0 eta1 eta2 eta3'
    sum_materials = 'fch_sn2 fe2 fp2'
    use_displaced_mesh = true
  [../]
  [./F3]
    type = DerivativeSumMaterial
    f_name = F3
    args = 'c3 eta0 eta1 eta2 eta3'
    sum_materials = 'fch_imc fe3'
    use_displaced_mesh = true
  [../]

  ## MOBILITIES
  # Diffusion coefficients
  [./diffusion_coeff]
    type = GenericConstantMaterial
    prop_names = 'Dcu Dimc Dsn'
    prop_values = '4.04e-26 8.11e-16 8.05e-15'
  [../]
  [./diffusion_coeff_gb]
    type = ParsedMaterial
    material_property_names = 'Dsn'
    f_name = Dgb
    function = '200*Dsn'
    use_displaced_mesh = true
  [../]

  #CahnHilliard Mobility
  [./Mbulk]
    type = DerivativeParsedMaterial
    f_name = M
    args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'length_scale energy_scale h0 h1 h2 h3 Dcu Dimc Dsn A_cu A_imc A_sn'
    function = '(length_scale^5/energy_scale)*(h0*Dcu/A_cu+(h1+h2)*Dsn/A_sn+h3*Dimc/A_imc)'
    use_displaced_mesh = true
  [../]
  # TODO GRAIN BOUNDARY MOBILITY NOT ADDED YET
  # AllenCahn mobility
  # [./ACMobility]
  #     type = GenericConstantMaterial
  #     prop_names = L
  #     prop_values = 2.7 #2.7
  #     use_displaced_mesh = true
  # [../]
  [./Lcuimc]
    type = ParsedMaterial
    f_name = Lcuimc
    material_property_names = 'energy_scale length_scale chat_cu chat_imc Dcu Dimc A_cu A_imc kappa mu'
    # function = '(length_scale^5/energy_scale)*mu*(Dcu/A_cu+Dimc/A_imc)/(3*kappa*(chat_cu-chat_imc)^2)'
    function = '2.7'
    output_properties = Lcuimc
    outputs = exodus
  [../]
  [./Lsnimc]
    type = ParsedMaterial
    f_name = Lsnimc
    material_property_names = 'energy_scale length_scale chat_sn chat_imc Dsn Dimc A_sn A_imc kappa mu'
    # function = '(length_scale^5/energy_scale)*mu*(Dsn/A_sn+Dimc/A_imc)/(3*kappa*(chat_sn-chat_imc)^2)'
    function = '2.7'
    output_properties = Lsnimc
    outputs = exodus
  [../]
  [./Lsnsn]
    type = ParsedMaterial
    f_name = Lsnsn
    material_property_names = 'energy_scale length_scale chat_sn chat_imc Dsn Dimc A_sn A_imc kappa mu'
    # function = '(length_scale^5/energy_scale)*mu*(Dsn/A_sn+Dimc/A_imc)/(3*kappa*(chat_sn-chat_imc)^2)'
    function = '2.7'
    output_properties = Lsnsn
    outputs = exodus
  [../]
  # [./ACMobility]
  #     type = DerivativeParsedMaterial
  #     args = 'eta0 eta1 eta2 eta3 eta4'
  #     f_name = L
  #     material_property_names = 'Lcuimc Lsnimc Lsnsn'
  #     constant_names = 'pf eps'
  #     constant_expressions = '1 0'
  #     # function='denom:=(pf*eta0^2+eps)*((pf*eta1^2+eps)+(pf*eta2^2+eps)+(pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta2^2+eps)*((pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta3^2+eps)*((pf*eta4^2+eps)); (Lcuimc*((pf*eta0^2+eps)*(pf*eta4^2+eps))+Lsnimc*((pf*eta4^2+eps))*((pf*eta1^2+eps)+(pf*eta2^2+eps)+(pf*eta3^2+eps))+Lsnsn*((pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+(pf*eta2^2+eps)*((pf*eta3^2+eps))))/denom'
  #     function='denom:=(pf*eta0^2+eps)*((pf*eta1^2+eps)+(pf*eta2^2+eps)+(pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta2^2+eps)*((pf*eta3^2+eps)+(pf*eta4^2+eps))+(pf*eta3^2+eps)*((pf*eta4^2+eps)); if(denom>1e-6,(Lcuimc*((pf*eta0^2+eps)*(pf*eta4^2+eps))+Lsnimc*((pf*eta4^2+eps))*((pf*eta1^2+eps)+(pf*eta2^2+eps)+(pf*eta3^2+eps))+Lsnsn*((pf*eta1^2+eps)*((pf*eta2^2+eps)+(pf*eta3^2+eps))+(pf*eta2^2+eps)*((pf*eta3^2+eps))))/denom,0)'
  #     # function='2.7'
  #     use_displaced_mesh = true
  #     output_properties = L
  #     outputs = exodus
  # [../]
  [./ACMobility]
    type = GenericConstantMaterial
    prop_names = L
    # prop_values = 2.7
    prop_values = 2.7
  [../]

  ## MODEL PARAMETERS
  [./kappa]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = kappa
    function = '0.75*sigma*delta*energy_scale/length_scale' #eV/nm
    use_displaced_mesh = true
  [../]
  [./mu]
    type = ParsedMaterial
    material_property_names = 'sigma delta length_scale energy_scale'
    f_name = mu
    function = '6*(sigma/delta)*energy_scale/length_scale^3' #eV/nm^3
    output_properties = mu
    use_displaced_mesh = true
  [../]
[]

[Kernels]
  #Stress divergence
  [./TensorMechanics]
    #displacements = 'disp_x disp_y disp_z'
    base_name = global
    use_displaced_mesh = true
    strain = FINITE #this also sets incremental strain =true
  [../]
  #Cahn-Hilliard equation
  [./CHC]
    type = KKSSplitCHCRes
    variable = c
    ca = c3
    cb = c1
    fa_name = F3
    fb_name = F1
    w = w
    h_name = h3
    args_a = 'eta0 eta1 eta2 eta3'
    use_displaced_mesh = true
  [../]
  [./dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
    use_displaced_mesh = true
  [../]
  [./CHW]
    type = SplitCHWRes
    mob_name = M # TODO CHANGE WHEN GB MOBILITY GETS ADDED
    variable = w
    args = 'eta0 eta1 eta2 eta3'
    use_displaced_mesh = true
  [../]
  # KKS conditions
  [./mu_cu_imc]
    type = KKSPhaseChemicalPotential
    variable = c0
    cb = c3
    fa_name = F0
    fb_name = F3
    args_a = ' '
    args_b = 'eta0 eta1 eta2 eta3'
    use_displaced_mesh = true
  [../]
  [./mu_imc_sn]
    type = KKSPhaseChemicalPotential
    variable = c3
    cb = c2
    fa_name = F3
    fb_name = F2
    args_a = 'eta0 eta1 eta2 eta3'
    args_b = 'eta0 eta1 eta2 eta3 '
    use_displaced_mesh = true
  [../]
  [./mu_sn_sn]
    type = KKSPhaseChemicalPotential
    variable = c2
    cb = c1
    fa_name = F2
    fb_name = F1
    args_a = 'eta0 eta1 eta2 eta3 '
    args_b = 'eta0 eta1 eta2 eta3 '
    use_displaced_mesh = true
  [../]
  [./mu_sn_sn2]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb = c0
    fa_name = F1
    fb_name = F0
    args_a = 'eta0 eta1 eta2 eta3 '
    args_b = ' eta0 eta1 eta2 eta3'
    use_displaced_mesh = true
  [../]
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = c0
    c = c
    cj = 'c0 c1 c2 c3'
    hj_names = 'h0 h1 h2 h3 '
    etas = 'eta0 eta1 eta2 eta3'
    use_displaced_mesh = true
  [../]
  #Allen-Cahn equations
  # AC Kernels
  [./KKSMultiACKernel]
    op_num = 4
    op_name_base = 'eta'
    ci_name_base = 'c'
    f_name_base = F
    #wi = 0.0624
    #wi = 10.
    #wi = 4.
    wi = 0.
    g_name_base = g
    use_displaced_mesh = true
  [../]

[]

[AuxVariables]
  [./imc_sn]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./imc_cu]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  [./shyd]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_tot]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slips]
    order = FIRST
    family = MONOMIAL
  [../]
  [./accum_slip]
    order = FIRST
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./slips]
    type =  CrystalPlasticityTotalSlip
    variable = slips
    h_names = 'h1 h2'
    slip_rates = 'slip_rate1 slip_rate2'
    execute_on = 'initial TIMESTEP_END'
  [../]
  [./accumslip]
    type = AccumulateAux
    accumulate_from_variable = slips
    variable = accum_slip
    execute_on = 'initial TIMESTEP_END'
  [../]
  [./f_tot]
    type = KKSMultiFreeEnergy
    variable = f_tot
    hj_names = 'h0 h1 h2 h3'
    Fj_names = 'F0 F1 F2 F3'
    gj_names = 'g0 g1 g2 g3'
    interfacial_vars = 'eta0 eta1 eta2 eta3'
    kappa_names = 'kappa kappa kappa kappa'
    additional_free_energy = f_int
    w = 0.
    use_displaced_mesh = true
  [../]
  [./f_int]
    type = ParsedAux
    variable = f_int
    args = 'eta0 eta1 eta2 eta3'
    constant_names = 'sigma delta gamma length_scale energy_scale'
    constant_expressions = '0.5 50e-9 1.5 1e9 6.24150943e18'
    function = 'mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta0^4-0.5*eta0^2+0.25*eta1^4-0.5*eta1^2+0.25*eta2^4-0.5*eta2^2+0.25*eta3^4-0.5*eta3^2+0.5*gamma*(eta0^2*(eta1^2+eta2^2+eta3^2)+eta1^2*(eta0^2+eta2^2+eta3^2)+eta2^2*(eta0^2+eta1^2+eta3^2)+eta3^2*(eta0^2+eta1^2+eta2^2)))'
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./imc_sn]
    type = ParsedAux
    variable = imc_sn
    args = 'h1 h2 h3 imc_sn'
    function = 'if(h3*(h1+h2)>imc_sn,h3*(h1+h2),imc_sn)'
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  [../]
  [./imc_cu]
    type = ParsedAux
    variable = imc_cu
    args = 'h0 h3 imc_cu'
    function = 'if(h3*h0>imc_cu,h3*h0,imc_cu)'
  [../]
  [./shyd]
    type = RankTwoScalarAux
    rank_two_tensor = global_stress
    variable = shyd
    scalar_type = Hydrostatic
  [../]
[]

[VectorPostprocessors]
  [./slip_line]
    type = LineValueSampler
    contains_complete_history = false
    variable = accum_slip
    start_point = '-500 0 0'
    end_point = '500 0 0'
    num_points = 201
    sort_by = x
    use_displaced_mesh = true
  [../]
  [./eta1_line]
    type = LineValueSampler
    contains_complete_history = false
    variable = eta1
    start_point = '-500 0 0'
    end_point = '500 0 0'
    num_points = 201
    sort_by = x
    use_displaced_mesh = true
  [../]
  [./eta2_line]
    type = LineValueSampler
    contains_complete_history = false
    variable = eta2
    start_point = '-500 0 0'
    end_point = '500 0 0'
    num_points = 201
    sort_by = x
    use_displaced_mesh = true
  [../]

[]
[Postprocessors]
  [./imc_area]
    type = ElementIntegralMaterialProperty
    mat_prop = h3
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_tot
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./step_size]
    type = TimestepSize
  [../]
[]


[Debug]
  show_var_residual_norms = true
  show_material_props = false
[]
[Preconditioning]
  [./smp]
    type = SMP
    full = true
    solve_type = PJFNK
    petsc_options_iname = '-pc_asm_overlap -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type -sub_ksp_type'
    petsc_options_value = '2                  asm       lu             nonzero             nonzero                    preonly'
  [../]
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options = '-snes_converged_reason -ksp_converged_reason
  -snes_ksp_ew'
  line_search = basic


  end_time = 100
  # num_steps = 2
  nl_max_its = 15
  l_max_its = 25

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_tol = 1e-3

  [./TimeIntegrator]
      type = LStableDirk2
  [../]
  [./TimeStepper]
      type = IterationAdaptiveDT
      dt =5e-3
      cutback_factor = 0.5
      growth_factor = 1.5
      optimal_iterations = 16 # 20
      linear_iteration_ratio = 25
  [../]
[]
[Outputs]
  # file_base = circle_100x100_0-0-0_90-90-0_0.003_symall_hex_fp_wi0
  file_base = circle2d_200x100_0.02_sym_fp_wi0_0-0-0_90-90-0_half
  [./exodus]
    type = Exodus
    append_date = true
  [../]
  [./csv]
    type = CSV
    append_date = true
  [../]

  perf_graph = true
[]

[Problem]
  solve = true
[]
