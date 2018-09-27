[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 80
  ny = 114
  nz = 1
  xmin = 0
  xmax = 1280
  ymin = -824
  ymax = 1000
  zmin = 0
  zmax = 16
  displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]
  # CahnHilliard needs the third derivatives
  derivative_order = 3
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    scaling = 1e-2
  [../]
  [./disp_y]
    scaling = 1e-2
  [../]
  [./disp_z]
    scaling = 1e-2
  [../]
  [./c]
    #scaling = 1e2
  [../]
  # chemical potential
  [./w]
    #scaling = 1e1
    #initial_condition = -0.7405
  [../]
  #Cu
  [./eta0]
    #scaling = 1e1
  [../]
  # Cu6Sn5
  [./eta1]
    initial_condition = 0
    #scaling = 1e1
  [../]
  #Sn
  [./eta2]
    #initial_condition = 1
    #scaling = 1e1
  [../]
  [./eta3]
    #initial_condition = 1
    #scaling = 1e1
  [../]
  #Cu
  [./c0]
    initial_condition = 0.02
    #initial_condition = 0.1617
    #scaling = 1e2
  [../]
  #Cu6Sn5
  [./c1]
    initial_condition = 0.4350
    #scaling = 1e2
  [../]
  #Sn
  [./c2]
    initial_condition = 0.9745
    #scaling = 1e2
  [../]
  [./c3]
    initial_condition = 0.9745
    #scaling = 1e2
  [../]
[]

[ICs]
  [./eta2] #Central Sn grain
    type = BoundingBoxIC
    variable = eta2
    inside = 1
    outside = 0
    x1 = 320 #100
    y1 = 0
    z1 = 0
    x2 = 960 #300
    y2 = 1000
    z2 = 16
    # type = FunctionIC
    # variable = eta2
    # function = Sn_grain_2_IC_fun
  [../]
  [./eta0] #Cu grain
    type = FunctionIC
    variable = eta0
    #function = 'if(y<1000,1,0)' #abaqus
    function = 'if(y<0,1,0)' #generated
  [../]
  [./eta3]
    type = UnitySubVarIC
    variable = eta3
    etas = 'eta2 eta0'
  [../]
  [./c] #Concentration of Sn
    type = VarDepIC
    variable = c
    cis = 'c0 c1 c2 c3'
    etas = 'eta0 eta1 eta2 eta3'
  [../]
  [./w]
    type = FunctionIC
    variable = w
    #function = 'if(y<1000,-16.4486,-0.74034)' #abaqus
    function = 'if(y<0,-16.4486,-0.74034)' #generated
  [../]
[]

[BCs]
  [./Periodic]
    #generated mesh
    [./xy]
      auto_direction = 'x z'
      variable = 'eta0 eta1 eta2 eta3 c w c0 c1 c2 c3 disp_x disp_y disp_z'
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
    boundary = 'bottom'# left right'
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = 'bottom'# front back'
    value = 0
  [../]
[]

[Functions]
  [./tanh_arg]
    type = ParsedFunction
    vars = 't_start_cool t_end_cool k_tan'
    vals = '1 21 0.333333'
    value = 'k_tan*(t - (t_start_cool + t_end_cool)/2)'
  [../]
  # [./Sn_grain_2_IC_fun] #sym line x_max/2 two lines cutting of the mid grain DEPENDS ON THE MESH
  #   type = ParsedFunction
  #   vars = 'x_11  x_12  y_11  y_12  m_1       m_2'
  #   vals = '310   320   1000  0     3.2e4     -9.6e4'
  #   value = 'if(y > 0 & (y > (((y_12-y_11)/(x_12-x_11)*x + m_1)) & y > ((-(y_12-y_11)/(x_12-x_11)*x + m_2))), 1, 0)'
  # [../]
[]

[UserObjects]
  #Crystal plasticity for central Sn grain
  [./slip_rate_gss2]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss2
    base_name = 'eta2'
  [../]
  [./slip_resistance_gss2]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss2
  [../]
  [./state_var_gss2]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    #groups = '0 32'
    #group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss2
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss2]
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

    uo_slip_rate_name = slip_rate_gss2
    uo_state_var_name = state_var_gss2
  [../]
  #Crystal plasticity for side Sn grain
  [./slip_rate_gss3]
    type = CrystalPlasticitySlipRateGSSBaseName
    variable_size = 32
    slip_sys_file_name = slip_systems_bct.txt
    num_slip_sys_flowrate_props = 2
    #flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1' #start_ss end_ss gamma0 1/m
    flowprops = '1 32 0.001 0.166666667' #start_ss end_ss gamma0 1/m
    uo_state_var_name = state_var_gss3
    base_name = 'eta3'
  [../]
  [./slip_resistance_gss3]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 32
    uo_state_var_name = state_var_gss3
  [../]
  [./state_var_gss3]
    type = CrystalPlasticityStateVariable
    variable_size = 32
    #groups = '0 32'
    #group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
    groups = '0 2 4 6 10 12 16 18 20 24 32'
    # initial slip resistances
    group_values = '0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss3
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss3]
    type = CrystalPlasticityStateVarRateComponentVoceP
    variable_size = 32
    #hprops = '1.4 100 40 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    #hprops = '1.4 0.624 0.250 2' #qab h0 ss c see eq (9) in Zhao 2017, values from Darbandi 2013 table V
    groups = '0 2 4 6 10 20 24 32'
    h0_group_values = '0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484'
    tau0_group_values = '0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    tauSat_group_values = '0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115'
    hardeningExponent_group_values = '2.0 2.0 2.0 2.0 2.0 2.0 2.0'
    coplanarHardening_group_values = '1.0 1.0 1.0 1.0 1.0 1.0 1.0' #q_aa = 1
    selfHardening_group_values = '1.4 1.4 1.4 1.4 1.4 1.4 1.4'

    crystal_lattice_type = BCT


    uo_slip_rate_name = slip_rate_gss3
    uo_state_var_name = state_var_gss3
  [../]
  # [./perf_log_dump]
  #   type = PerflogDumper
  #   execute_on = TIMESTEP_END
  #   outfile = perflog.csv
  #   use_displaced_mesh = true
  # [../]
[]

[Materials]
  [./time]
    type = TimeStepMaterial
    prop_time = time
    prop_dt = dt
    use_displaced_mesh = true
    # output_properties = 'dt'
    # outputs = exodus
  [../]
  # Cu
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
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    base_name = eta0
    f_name = fe0
    outputs = exodus
    output_properties = fe0
    use_displaced_mesh = true
  [../]
  #central Sn
  [./crysp2]
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate_gss2'
    uo_slip_resistances = 'slip_resistance_gss2'
    uo_state_vars = 'state_var_gss2'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss2'
    base_name = 'eta2'
    maximum_substep_iteration = 10
    tan_mod_type = exact
    # output_properties = 'slip_rate_gss2'
    # outputs = exodus
  [../]
  [./strain2]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = 'eta2'
  [../]
  [./elasticity_tensor2]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3' #MPa #From Darbandi 2014 table I
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    base_name = 'eta2'
  [../]
  [./fe2]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe2
    base_name = eta2
    use_displaced_mesh = true #Not sure
    plasticity = true
    outputs = exodus
    output_properties = fe2
  [../]
  [./fp2]
    type = CPPlasticEnergyMaterial
    # Q = 1 in class
    variable_size = 32
    uo_state_var_name = state_var_gss2
    f_name = fp2
    use_displaced_mesh = true
    args = ' '
    #s0 = 0.144
    groups = '0 2 4 6 10 12 16 18 20 24 32' #values calibrated on 110
    group_values = '0.05122774 0.03452174 0.03696807 0.00912421 0.02046358 0.01612225 0.04525029 0.08612754 0.29181706 0.02277457'
    # outputs = exodus
    # output_properties = fp2
  [../]
  #Side Sn
  [./crysp3]
    type = FiniteStrainUObasedCPBaseName
    rtol = 1e-6
    abs_tol = 1e-6
    stol = 1e-2
    uo_slip_rates = 'slip_rate_gss3'
    uo_slip_resistances = 'slip_resistance_gss3'
    uo_state_vars = 'state_var_gss3'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss3'
    base_name = 'eta3'
    maximum_substep_iteration = 10
    tan_mod_type = exact
    # output_properties = 'slip_rate_gss3'
    # outputs = exodus
  [../]
  [./strain3]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = 'eta3'
  [../]
  [./elasticity_tensor3]
    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation
    #C_ijkl = '72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 24e3 22e3 22e3' #MPa #From Darbandi 2013 table III
    C_ijkl = '451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022' #eV/nm^3 #From Darbandi 2013 table I
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    base_name = 'eta3'
  [../]
  [./fe3]
    type = ElasticEnergyMaterialGreenPK2
    args = ' '
    f_name = fe3
    base_name = eta3
    use_displaced_mesh = true #Not sure
    plasticity = true
    outputs = exodus
    output_properties = fe3
  [../]
  [./fp3]
    type = CPPlasticEnergyMaterial
    # Q = 1 in class
    variable_size = 32
    uo_state_var_name = state_var_gss3
    f_name = fp3
    use_displaced_mesh = true
    args = ' '
    #s0 = 0.144
    groups = '0 2 4 6 10 12 16 18 20 24 32' #values calibrated on 110
    group_values = '0.05122774 0.03452174 0.03696807 0.00912421 0.02046358 0.01612225 0.04525029 0.08612754 0.29181706 0.02277457'
    # outputs = exodus
    # output_properties = fp3
  [../]

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
    args = 'eta0 eta1 eta2 eta3'
    #args = ' '
    base_name = eta1
    eigen_base = '1 1 1 0 0 0'
    #eigen_base = '0.02 0.02 0'
    #eigen_base = '0 0 0'
    eigenstrain_name = eT_eta
    prefactor = pre
  [../]
  [./pre]
    type = DerivativeParsedMaterial
    args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'h1'
    function = '-0.003*h1'
    f_name = pre
    outputs = exodus
  [../]
  [./fe1]
    type = ElasticEnergyMaterialGreenPK2
    args = 'eta0 eta1 eta2 eta3'
    #args = ' '
    f_name = fe1
    base_name = eta1
    eigenstrain = true
    eigenstrain_name = eT_eta
    use_displaced_mesh = true #Not sure
    outputs = exodus
    output_properties = fe1
  [../]
  [./scale]
    type = GenericConstantMaterial
    prop_names = 'length_scale energy_scale time_scale'
    prop_values = '1e9 6.24150943e18 1.' #m to nm J to eV s to h
  [../]
  [./model_constants]
    type = GenericConstantMaterial # delta is diff reg. width
    prop_names = 'sigma delta delta_real gamma tgrad_corr_mult'
    prop_values = '0.5 80e-9 5e-10 1.5 0' #J/m^2 m - ?
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
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h0
  [../]
  [./h1]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h1
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta1
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h1
  [../]
  [./h2]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h2
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta2
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h2
  [../]
  [./h3]
      type = SwitchingFunctionMultiPhaseMaterial
      h_name = h3
      all_etas = 'eta0 eta1 eta2 eta3'
      phase_etas = eta3
      use_displaced_mesh = true
      outputs = exodus
      output_properties = h3
  [../]
  [./g0]
    type = BarrierFunctionMaterial
    eta = eta0
    well_only = true
    function_name = g0
    g_order = SIMPLE
    use_displaced_mesh = true
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
  [./g3]
    type = BarrierFunctionMaterial
    eta = eta3
    well_only = true
    function_name = g3
    g_order = SIMPLE
    use_displaced_mesh = true
  [../]
  [./ACMobility]
      type = GenericConstantMaterial
      prop_names = L
      prop_values = 2.7 #2.7
  [../]
  [./noise_constants]
    type = GenericConstantMaterial
    prop_names = 'T kb lambda dim' #eta2erature Boltzmann gridsize dimensionality
    #prop_values = '493 8.6173303e-5  3'
    prop_values = '493 8.6173303e-5 8 3'
  [../]
  [./nuc]
    type =  DerivativeParsedMaterial
    args = 'eta0 eta1 eta2 eta3'
    f_name = nuc
    material_property_names = 'time dt T kb lambda dim L h0 h2 h3'
    function = 'if(time<1&(h0*h2>0.09 | h0*h3>0.09),sqrt(2*kb*T*L/(lambda^dim*dt)),0)' #expression from Shen (2007) not sure about diff reg coeff hi*hj>X?
    outputs = exodus
    output_properties = nuc
    use_displaced_mesh = true
  [../]
  # Constants The energy parameters are for 220 C and for 25 C =================
  [./energy_constants_A]
    type =  GenericFunctionMaterial # GenericConstantMaterial
    prop_names = 'A_cu_220  A_eps_220   A_eta_220   A_sn_220
                  A_cu_25   A_eps_25    A_eta_25    A_sn_25'
    #prop_values = '1.0133e5/Vm 4e5/Vm 4.2059e6/Vm' #J/m^3
    prop_values = '1.7756e10  2.4555e11   2.4555e11   2.3033e10
                   6.2204e9   2.4555e10   2.4555e10   2.5819e11' #J/m^3 Aeps=Aeta=2e6 prev. used; Vm := molar vol = 16.29 [cm^3 mol^-1]
    #prop_values = '1.5929e10 2.4555e12 2.4555e12 2.3020e10' #J/m^3 Aeps = 2e7 Aeta = 2e7
  [../]
  [./energy_constants_B]
    type = GenericConstantMaterial
    prop_names = 'B_cu_220  B_eps_220 B_eta_220 B_sn_220
                  B_cu_25   B_eps_25  B_eta_25  B_sn_25'
    #prop_values = '-2.1146e4/Vm -6.9892e3/Vm 7.168e3/Vm' #J/m^3
    prop_values = '-2.6351e9  -1.4014e9   2.3251e7    2.14216e8
                   -1.2981e9  -4.2905e8   -4.2905e8   4.4002e8' #J/m^3
    #prop_values = '-2.5789e9 -1.3733e9 2.3175e7 2.1406e8' #J/m^3
  [../]
  [./energy_C]
    type = GenericConstantMaterial
    prop_names = 'C_cu_220  C_eps_220  C_eta_220  C_sn_220
                  C_cu_25   C_eps_25   C_eta_25   C_sn_25'
    #prop_values = '-1.2842e4/Vm -1.9185e4/Vm -1.5265e4/Vm' #J/m^3
    prop_values = '-1.1441e9  -1.7294e9   -1.7646e9   -1.646e9
                   -7.8834e8  -1.1778e9   -1.1778e9   -9.3708e8' #J/m^3
    #prop_values = '-1.1529e9 -1.7330e9 -1.7646e9 -1.646e9' #J/m^3
  [../]
  # [./sim_time_constants]
  #   type = GenericConstantMaterial
  #   #times for the ramping of mat. constants
  #   prop_names = 't_start_cool t_end_cool'
  #   prop_values = '3 23' # [s]
  # [../]

  [./tanh_arg]
    type = GenericFunctionMaterial
    #times for the ramping of mat. constants
    prop_names = 'phi_t'
    prop_values = 'tanh_arg' # [s]
  [../]
  #====================== GIBBS ENERGY CONSTANTS ABC ===========================
  [./energy_A_cu]
      type =  ParsedMaterial
      material_property_names = 'A_cu_220 A_cu_25 phi_t'
      function = 'A_cu_25 - (A_cu_25 - A_cu_220)*(1-tanh(phi_t))/2'
      f_name = A_cu
  [../]
  [./energy_A_eps]
      type =  ParsedMaterial
      material_property_names = 'A_eps_220 A_eps_25 phi_t'
      function = 'A_eps_25 - (A_eps_25 - A_eps_220)*(1-tanh(phi_t))/2'
      f_name = A_eps
  [../]
  [./energy_A_eta]
      type =  ParsedMaterial
      material_property_names = 'A_eta_220 A_eta_25 phi_t'
      function = 'A_eta_25 - (A_eta_25 - A_eta_220)*(1-tanh(phi_t))/2'
      f_name = A_eta
  [../]
  [./energy_A_sn]
      type =  ParsedMaterial
      material_property_names = 'A_sn_220 A_sn_25 phi_t'
      function = 'A_sn_25 - (A_sn_25 - A_sn_220)*(1-tanh(phi_t))/2'
      f_name = A_sn
  [../]
  [./energy_B_cu]
      type =  ParsedMaterial
      material_property_names = 'B_cu_220 B_cu_25 phi_t'
      function = 'B_cu_25 - (B_cu_25 - B_cu_220)*(1-tanh(phi_t))/2'
      f_name = B_cu
  [../]
  [./energy_B_eps]
      type =  ParsedMaterial
      material_property_names = 'B_eps_220 B_eps_25 phi_t'
      function = 'B_eps_25 - (B_eps_25 - B_eps_220)*(1-tanh(phi_t))/2'
      f_name = B_eps
  [../]
  [./energy_B_eta]
      type =  ParsedMaterial
      material_property_names = 'B_eta_220 B_eta_25 phi_t'
      function = 'B_eta_25 - (B_eta_25 - B_eta_220)*(1-tanh(phi_t))/2'
      f_name = B_eta
  [../]
  [./energy_B_sn]
      type =  ParsedMaterial
      material_property_names = 'B_sn_220 B_sn_25 phi_t'
      function = 'B_sn_25 - (B_sn_25 - B_sn_220)*(1-tanh(phi_t))/2'
      f_name = B_sn
  [../]
  [./energy_C_cu]
      type =  ParsedMaterial
      material_property_names = 'C_cu_220 C_cu_25 phi_t'
      function = 'C_cu_25 - (C_cu_25 - C_cu_220)*(1-tanh(phi_t))/2'
      f_name = C_cu
  [../]
  [./energy_C_eps]
      type =  ParsedMaterial
      material_property_names = 'C_eps_220 C_eps_25 phi_t'
      function = 'C_eps_25 - (C_eps_25 - C_eps_220)*(1-tanh(phi_t))/2'
      f_name = C_eps
  [../]
  [./energy_C_eta]
      type =  ParsedMaterial
      material_property_names = 'C_eta_220 C_eta_25 phi_t'
      function = 'C_eta_25 - (C_eta_25 - C_eta_220)*(1-tanh(phi_t))/2'
      f_name = C_eta
  [../]
  [./energy_C_sn]
      type =  ParsedMaterial
      material_property_names = 'C_sn_220 C_sn_25 phi_t'
      function = 'C_sn_25 - (C_sn_25 - C_sn_220)*(1-tanh(phi_t))/2'
      f_name = C_sn
  [../]

  #=================== END OF GIBBS ENERGY PARAMETER CONSTANTS ABC =============
  # [./energy_c_ab]
  #   type = GenericConstantMaterial
  #   prop_names = 'c_cu_eps c_cu_eta c_cu_sn c_eps_cu c_eps_eta c_eps_sn c_eta_cu c_eta_eps c_eta_sn c_sn_cu c_sn_eps c_sn_eta'
  #   prop_values = '0.02 0.1957 0.6088 0.2383 0.2483 0.2495 0.4299 0.4343 0.4359 0.9789 0.9839 0.9889' #-
  #   #prop_values = '0.0234 0.198 0.6088 0.2479 0.2489 0.000 0.4345 0.4349 0.4351 0.9789 0.000 0.9889' #- Aeps = 2e7 Aeta = 2e7
  # [../]
  #============ GIBBS ENERGY EQ. MOLAR CONC. =================================== eps phase is redundant
  [./energy_constants_chat]
    type = GenericConstantMaterial
    prop_names = 'chat_cu_220   chat_eps_220  chat_eta_220  chat_sn_220
                  chat_cu_25    chat_eps_25   chat_eta_25   chat_sn_25'
    prop_values = '0.02     0.2433    0.4351    0.9889
                   0.10569  0.41753   0.41753   0.99941' #-
    #prop_values = '0.0234 0.2484 0.4350 0.9889' #-
  [../]
  [./diffusion_constants]
    type = GenericConstantMaterial
    prop_names = 'D_cu_220  D_eps_220   D_eta_220   D_sn_220
                  D_cu_25   D_eps_25    D_eta_25    D_sn_25'
    #prop_values = '1e-20 6e-16 1.5e-14 1e-13' # m^2/s #D12 best slightly slow
    #prop_values = '1e-20 9.5e-16 3e-14 1e-13' # m^2/s #D15
    prop_values = '1e-20      1.25e-15  3.1e-14   1e-13
                   2.877e-36  6.575e-19 6.575e-19 2.452e-17' # m^2/s #D16 BEST
    #prop_values = '1e-16 1.25e-15 3.1e-14 1e-13' # m^2/s #D17
    #outputs = exodus
  [../]
  [./energy_chat_cu]
      type =  ParsedMaterial
      material_property_names = 'chat_cu_220 chat_cu_25 phi_t'
      function = 'chat_cu_25 - (chat_cu_25 - chat_cu_220)*(1-tanh(phi_t))/2'
      f_name = chat_cu
  [../]
  [./energy_chat_eps]
      type =  ParsedMaterial
      material_property_names = 'chat_eps_220 chat_eps_25 phi_t'
      function = 'chat_eps_25 - (chat_eps_25 - chat_eps_220)*(1-tanh(phi_t))/2'
      f_name = chat_eps
  [../]
  [./energy_chat_eta]
      type =  ParsedMaterial
      material_property_names = 'chat_eta_220 chat_eta_25 phi_t'
      function = 'chat_eta_25 - (chat_eta_25 - chat_eta_220)*(1-tanh(phi_t))/2'
      f_name = chat_eta
  [../]
  [./energy_chat_sn]
      type =  ParsedMaterial
      material_property_names = 'chat_sn_220 chat_sn_25 phi_t'
      function = 'chat_sn_25 - (chat_sn_25 - chat_sn_220)*(1-tanh(phi_t))/2'
      f_name = chat_sn
  [../]
  [./bulk_diff_D_cu]
      type =  ParsedMaterial
      material_property_names = 'D_cu_220 D_cu_25 phi_t'
      function = 'D_cu_25 - (D_cu_25 - D_cu_220)*(1-tanh(phi_t))/2'
      f_name = D_cu
  [../]
  [./bulk_diff_D_eps]
      type =  ParsedMaterial
      material_property_names = 'D_eps_220 D_eps_25 phi_t'
      function = 'D_eps_25 - (D_eps_25 - D_eps_220)*(1-tanh(phi_t))/2'
      f_name = D_eps
  [../]
  [./bulk_diff_D_eta]
      type =  ParsedMaterial
      material_property_names = 'D_eta_220 D_eta_25 phi_t'
      function = 'D_eta_25 - (D_eta_25 - D_eta_220)*(1-tanh(phi_t))/2'
      f_name = D_eta
  [../]
  [./bulk_diff_D_sn]
      type =  ParsedMaterial
      material_property_names = 'D_sn_220 D_sn_25 phi_t'
      function = 'D_sn_25 - (D_sn_25 - D_sn_220)*(1-tanh(phi_t))/2'
      f_name = D_sn
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
      use_displaced_mesh = true
      outputs = exodus
      output_properties = fch0
  [../]
  [./fch_imc] #Chemical energy Cu6Sn5
      type = DerivativeParsedMaterial
      f_name = fch1
      args = 'c1'
      material_property_names = 'A_eta B_eta C_eta chat_eta length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_eta*(c1-chat_eta)^2+B_eta*(c1-chat_eta)+C_eta)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
      outputs = exodus
      output_properties = fch1
  [../]
  [./fch_sn] #Chemical energy Sn central grain
      type = DerivativeParsedMaterial
      f_name = fch2
      args = 'c2'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c2-chat_sn)^2+B_sn*(c2-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
      outputs = exodus
      output_properties = fch2
  [../]
  [./fch_sn2] #Chemical energy Sn side grain
      type = DerivativeParsedMaterial
      f_name = fch3
      args = 'c3'
      material_property_names = 'A_sn B_sn C_sn chat_sn length_scale energy_scale'
      function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c3-chat_sn)^2+B_sn*(c3-chat_sn)+C_sn)' #eV/nm^3
      derivative_order = 2
      use_displaced_mesh = true
      outputs = exodus
      output_properties = fch3
  [../]

  [./Mgb]
    type=DerivativeParsedMaterial
    args = 'eta0 eta1 eta2 eta3'
    material_property_names = 'D_gb delta delta_real h0 h1 h2 h3 A_cu A_eta A_sn length_scale energy_scale time_scale'  # said A_imc instead of A_eta but no occurrence in code
    f_name = Mgb
    #function = '(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h1*A_imc+(h2+h3)*A_sn)*delta)'
    function = 'if(h2*h3>0.09,(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+h1*A_eta+(h2+h3)*A_sn)*delta),0)'
    #function = '4e-5'
    outputs = exodus
    output_properties = Mgb
  [../]
  [./Mbulk]
      type = DerivativeParsedMaterial
      f_name = Mbulk
      args = 'eta0 eta1 eta2 eta3'
      material_property_names = 'h0 h1 h2 h3 D_cu D_eta D_sn A_cu A_eta A_sn Mgb length_scale energy_scale time_scale'
      #function = 's:=eta_cu^2+eta_imc2^2+eta_imc1^2+eta_sn^2;p:=eta_imc2^2*eta_imc1^2;(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc2*D_imc/A_imc+h_imc1*D_imc/A_imc+h_sn*D_sn/A_sn+p*Mgb/s)' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc2*D_imc/A_imc+h_imc1*D_imc/A_imc+h_sn*D_sn/A_sn)+if(h_imc2*h_imc1>1./16.,0,Mgb)' #nm^5/eVs
      function = '(length_scale^5/(energy_scale*time_scale))*(h0*D_cu/A_cu+h1*D_eta/A_eta+(h2+h3)*D_sn/A_sn)' #nm^5/eVs
      #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
      derivative_order = 2
      #outputs = exodus
      #output_properties = M
      use_displaced_mesh = true
  [../]
  [./CHMobility]
    type = DerivativeSumMaterial
    f_name = M
    sum_materials = 'Mgb Mbulk'
    args = 'eta0 eta1 eta2 eta3'
    outputs = exodus
    output_properties = M
    use_displaced_mesh = true
  [../]
  [./F_cu]
    type = DerivativeSumMaterial
    f_name = F0
    args = 'c0'
    sum_materials = 'fch0 fe0'
    #sum_materials = 'fch2'
    use_displaced_mesh = true
  [../]
  [./F_eta]
    type = DerivativeSumMaterial
    f_name = F1
    args = 'c1 eta0 eta1 eta2 eta3'
    sum_materials = 'fch1 fe1'
    #sum_materials = 'fch1'
    use_displaced_mesh = true
  [../]
  [./F_sn2]
    type = DerivativeSumMaterial
    f_name = F2
    args = 'c2'
    sum_materials = 'fch2 fe2 fp2'
    #sum_materials = 'fch2'
    use_displaced_mesh = true
  [../]
  [./F_sn3]
    type = DerivativeSumMaterial
    f_name = F3
    args = 'c3'
    sum_materials = 'fch3 fe3 fp3'
    #sum_materials = 'fch2'
    use_displaced_mesh = true
  [../]
  # Generate the global stress from the phase stresses
  [./global_stress] #homogeniserar bara Cauchy stress
    type = MultiPhaseStressMaterial
    phase_base = 'eta0 eta1 eta2 eta3'
    h          = 'h0 h1 h2 h3'
    base_name = global
  [../]
  [./global_strain]
    type = ComputeFiniteStrain
    #displacements = 'disp_x disp_y disp_z'
    base_name = global
  [../]
[]

[Kernels]
  #Stress divergence
  [./TensorMechanics]
    #displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    base_name = global
    strain = FINITE #this also sets incremental strain =true
  [../]

  #Nucleation of Cu6Sn5
  [./nuceta1]
    type = LangevinNoisePositive
    variable = eta1
    amplitude = 1
    seed = 1e9 #123456789
    multiplier = nuc
    use_displaced_mesh = true
  [../]
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
      args_a = 'eta0 eta1 eta2 eta3'
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
      args = 'eta0 eta1 eta2 eta3'
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
    args_a = ' '
    args_b = 'eta0 eta1 eta2 eta3'
  [../]
  [./chempot_eta_sn]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb       = c2
    fa_name  = F1
    fb_name  = F2
    use_displaced_mesh = true
    args_a = 'eta0 eta1 eta2 eta3'
  [../]
  [./chempot_sn_sn]
    type = KKSPhaseChemicalPotential
    variable = c2
    cb       = c3
    fa_name  = F2
    fb_name  = F3
    use_displaced_mesh = true
    args_a = ' '
  [../]
  [./phaseconcentration] # enforce c = sum h_i*c_i
    type = KKSMultiPhaseConcentration
    variable = c3
    cj = 'c0 c1 c2 c3'
    hj_names = 'h0 h1 h2 h3'
    etas = 'eta0 eta1 eta2 eta3'
    c = c
    use_displaced_mesh = true
  [../]

  # AC Kernels
  [./KKSMultiACKernel]
    op_num = 4
    op_name_base = 'eta'
    ci_name_base = 'c'
    f_name_base = F
    #wi = 0.0624
    wi = 10.
    #wi = 4.
    #wi = 1.
    g_name_base = g
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
  [./hyd_sn3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hyd_sn3_mpa]
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
  [./f_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./hfp2]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./hfp3]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./sxx2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./szz2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxx3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./szz3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sbiax]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./step_dt]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./accum_L_2_slip_rates_gss_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./L_2_slip_rates_gss_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss21]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss23]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss24]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss25]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss26]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss27]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss28]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss29]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss210]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss211]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss212]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss213]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss214]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss215]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss216]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss217]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss218]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss219]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss220]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss221]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss222]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss223]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss224]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss225]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss226]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss227]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss228]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss229]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss230]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss231]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss232]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./accum_L_2_slip_rates_gss_3] #\sum d\gamma
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./L_2_slip_rates_gss_3] #d\gamma
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss31]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss32]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss33]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss34]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss35]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss36]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss37]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss38]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss39]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss310]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss311]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss312]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss313]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss314]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss315]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss316]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss317]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss318]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss319]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss320]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss321]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss322]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss323]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss324]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss325]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss326]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss327]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss328]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss329]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss330]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss331]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_rate_gss332]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./step_dt]
    type = MaterialRealAux
    property = dt
    variable = step_dt
    execute_on = TIMESTEP_END
  [../]
  # [./hfp2]
  #   type = ParsedAux
  #   variable = hfp2
  #   args = 'fp2 h2'
  #   function = 'h2*fp2'
  #   use_displaced_mesh = true
  # [../]
  # [./hfp3]
  #   type = ParsedAux
  #   variable = hfp3
  #   args = 'fp3 h3'
  #   function = 'h3*fp3'
  #   use_displaced_mesh = true
  # [../]
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
    rank_two_tensor = eta2_stress
    scalar_type = Hydrostatic
    execute_on = timestep_end
  [../]
  [./hyd_sn_mpa]
    type = ParsedAux
    variable = hyd_sn_mpa
    args = 'h2 hyd_sn'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h2*hyd_sn*to_MPa'
    execute_on = timestep_end
  [../]
  [./hyd_sn3]
    type = RankTwoScalarAux
    variable = hyd_sn3
    rank_two_tensor = eta1_stress
    scalar_type = Hydrostatic
    execute_on = timestep_end
  [../]
  [./hyd_sn3_mpa]
    type = ParsedAux
    variable = hyd_sn3_mpa
    args = 'h3 hyd_sn3'
    constant_names = 'to_MPa' #from eV/nm^3
    constant_expressions = '160.217662'
    function = 'h3*hyd_sn3*to_MPa'
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

  [./sxx2]
    type = RankTwoAux
    variable = sxx2
    rank_two_tensor = eta2_stress
    index_i = 0
    index_j = 0
  [../]
  [./szz2]
    type = RankTwoAux
    variable = szz2
    rank_two_tensor = eta2_stress
    index_i = 2
    index_j = 2
  [../]
  [./sxx3]
    type = RankTwoAux
    variable = sxx3
    rank_two_tensor = eta3_stress
    index_i = 0
    index_j = 0
  [../]
  [./szz3]
    type = RankTwoAux
    variable = szz3
    rank_two_tensor = eta3_stress
    index_i = 2
    index_j = 2
  [../]
  [./sbiax]
    type = ParsedAux
    variable = sbiax
    args = 'h2 h3 sxx2 szz2 sxx3 szz3'
    constant_names = 'to_MPa'
    constant_expressions = '160.217662'
    function = '(h2*0.5*(sxx2+szz2)+h3*0.5*(sxx3+szz3))*to_MPa'
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
    #w = 4.
    #w = 1.
    w = 0.
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./f_int]
    type = ParsedAux
    variable = f_int
    args = 'eta0 eta1 eta2 eta3'
    constant_names = 'sigma delta gamma length_scale energy_scale'
    constant_expressions = '0.5 80e-9 1.5 1e9 6.24150943e18'
    #constant_expressions = '0.5 60e-9 1.5 1e9 6.24150943e18'
    function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta0^4-0.5*eta0^2+0.25*eta2^4-0.5*eta2^2+0.25*eta1^4-0.5*eta1^2+gamma*(eta0^2*(eta2^2+eta1^2)+eta2^2*eta1^2)+0.25)'
    execute_on = 'initial timestep_end'
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss21]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss21
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss22]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss22
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss23]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss23
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss24]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss24
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss25]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss25
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss26]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss26
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss27]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss27
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss28]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss28
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss29]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss29
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss210]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss210
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss211]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss211
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss212]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss212
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss213]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss213
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss214]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss214
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss215]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss215
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss216]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss216
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss217]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss217
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss218]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss218
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss219]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss219
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss220]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss220
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss221]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss221
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss222]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss222
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss223]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss223
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss224]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss224
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss225]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss225
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss226]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss226
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss227]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss227
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss228]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss228
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss229]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss229
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss230]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss230
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss231]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss231
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss232]
    type = MaterialStdVectorAux
    property = slip_rate_gss2
    variable = slip_rate_gss232
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./L_2_slip_rates_gss_2]
    type = ParsedAux
      variable = L_2_slip_rates_gss_2
      args = 'slip_rate_gss21   slip_rate_gss22   slip_rate_gss23   slip_rate_gss24   slip_rate_gss25   slip_rate_gss26   slip_rate_gss27   slip_rate_gss28
              slip_rate_gss29   slip_rate_gss210   slip_rate_gss211   slip_rate_gss212   slip_rate_gss213   slip_rate_gss214   slip_rate_gss215   slip_rate_gss216
              slip_rate_gss217   slip_rate_gss218   slip_rate_gss219   slip_rate_gss220   slip_rate_gss221   slip_rate_gss222   slip_rate_gss223   slip_rate_gss224
              slip_rate_gss225   slip_rate_gss226   slip_rate_gss227   slip_rate_gss228   slip_rate_gss229   slip_rate_gss230   slip_rate_gss231   slip_rate_gss232
              h2 step_dt'
      function = 'step_dt*h2*sqrt(slip_rate_gss21^2 + slip_rate_gss22^2 +slip_rate_gss23^2 + slip_rate_gss24^2 + slip_rate_gss25^2 + slip_rate_gss26^2 + slip_rate_gss27^2 + slip_rate_gss28^2
                        + slip_rate_gss29^2 + slip_rate_gss210^2 + slip_rate_gss211^2 + slip_rate_gss212^2 + slip_rate_gss213^2 + slip_rate_gss214^2 + slip_rate_gss215^2 + slip_rate_gss216^2
                        + slip_rate_gss217^2 + slip_rate_gss218^2 + slip_rate_gss219^2 + slip_rate_gss220^2 + slip_rate_gss221^2 + slip_rate_gss222^2 + slip_rate_gss223^2 + slip_rate_gss224^2
                        + slip_rate_gss225^2 + slip_rate_gss226^2 + slip_rate_gss227^2 + slip_rate_gss228^2 + slip_rate_gss229^2 + slip_rate_gss230^2 + slip_rate_gss231^2 + slip_rate_gss232^2)'
  [../]
  [./accum_L_2_slip_rates_gss_2]
    type = AccumulateAux
    accumulate_from_variable = L_2_slip_rates_gss_2
    variable = accum_L_2_slip_rates_gss_2
    execute_on = TIMESTEP_END
  [../]
  [./slip_rate_gss31]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss31
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss32]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss32
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss33]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss33
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss34]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss34
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss35]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss35
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss36]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss36
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss37]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss37
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss38]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss38
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss39]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss39
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss310]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss310
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss311]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss311
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss312]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss312
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss313]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss313
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss314]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss314
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss315]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss315
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss316]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss316
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss317]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss317
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss318]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss318
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss319]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss319
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss320]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss320
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss321]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss321
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss322]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss322
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss323]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss323
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss324]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss324
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss325]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss325
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss326]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss326
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss327]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss327
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss328]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss328
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss329]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss329
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss330]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss330
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss331]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss331
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./slip_rate_gss332]
    type = MaterialStdVectorAux
    property = slip_rate_gss3
    variable = slip_rate_gss332
    index = 0
    execute_on = TIMESTEP_END
    use_displaced_mesh = true
  [../]
  [./L_2_slip_rates_gss_3]
    type = ParsedAux
      variable = L_2_slip_rates_gss_3
      args = 'slip_rate_gss31   slip_rate_gss32   slip_rate_gss33   slip_rate_gss34   slip_rate_gss35   slip_rate_gss36   slip_rate_gss37   slip_rate_gss38
              slip_rate_gss39   slip_rate_gss310   slip_rate_gss311   slip_rate_gss312   slip_rate_gss313   slip_rate_gss314   slip_rate_gss315   slip_rate_gss316
              slip_rate_gss317   slip_rate_gss318   slip_rate_gss319   slip_rate_gss320   slip_rate_gss321   slip_rate_gss322   slip_rate_gss323   slip_rate_gss324
              slip_rate_gss325   slip_rate_gss326   slip_rate_gss327   slip_rate_gss328   slip_rate_gss329   slip_rate_gss330   slip_rate_gss331   slip_rate_gss332
              h3 step_dt'
      function = 'step_dt*h3*sqrt(slip_rate_gss31^2 + slip_rate_gss32^2 +slip_rate_gss33^2 + slip_rate_gss34^2 + slip_rate_gss35^2 + slip_rate_gss36^2 + slip_rate_gss37^2 + slip_rate_gss38^2
                        + slip_rate_gss39^2 + slip_rate_gss310^2 + slip_rate_gss311^2 + slip_rate_gss312^2 + slip_rate_gss313^2 + slip_rate_gss314^2 + slip_rate_gss315^2 + slip_rate_gss316^2
                        + slip_rate_gss317^2 + slip_rate_gss318^2 + slip_rate_gss319^2 + slip_rate_gss320^2 + slip_rate_gss321^2 + slip_rate_gss322^2 + slip_rate_gss323^2 + slip_rate_gss324^2
                        + slip_rate_gss325^2 + slip_rate_gss326^2 + slip_rate_gss327^2 + slip_rate_gss328^2 + slip_rate_gss329^2 + slip_rate_gss330^2 + slip_rate_gss331^2 + slip_rate_gss332^2)'
  [../]
  [./accum_L_2_slip_rates_gss_3]
    type = AccumulateAux
    accumulate_from_variable = L_2_slip_rates_gss_3
    variable = accum_L_2_slip_rates_gss_3
    execute_on = TIMESTEP_END
  [../]
[]

[Postprocessors]
  [./imc_area_h]
    type = ElementIntegralMaterialProperty
    mat_prop = h1
    execute_on = 'INITIAL TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  [./sbiax_mean]
    type = ElementAverageValue # Should only calculate mean in Sn elements. Write a PhaseAverageValue maybe.
    variable = sbiax
    use_displaced_mesh = true
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = f_density
    execute_on = 'Initial TIMESTEP_END'
    use_displaced_mesh = true
  [../]
  [./step_size]
    type = TimestepSize
  [../]
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
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

   line_search = basic
   #line_search = none
   #line_search = default
   #line_search = bt

   petsc_options = '-snes_converged_reason -ksp_converged_reason
   -snes_ksp_ew'

   #l_max_its = 120
   l_max_its = 10
   #nl_max_its = 25
   nl_max_its = 15

   l_tol = 1.0e-3
   #l_abs_step_tol = 1e-8
   nl_rel_tol = 1.0e-8 # 1.0e-8 # 1.0e-9 #1.0e-10
   nl_abs_tol = 1.0e-9 # 1.0e-9 # 1.0e-10#1.0e-11

   #num_steps = 2000
   end_time = 1.2672e6 #2x176h
   #n_startup_steps = 2
   dtmin= 1e-5
   dtmax= 1e5

   #scheme = implicit-euler
   #scheme = bdf2
   # scheme = dirk
   [./TimeIntegrator]
     # type = ExplicitTVDRK2
       type = LStableDirk2
     # type = ImplicitEuler
   [../]

   [./TimeStepper]
       # Turn on time stepping
       type = IterationAdaptiveDT
       #dt = 1e-3
       dt = 5e-3


       cutback_factor = 0.5
       growth_factor = 1.3
       #optimal_iterations = 10
       optimal_iterations = 15 # 20
       linear_iteration_ratio = 25
       #postprocessor_dtlim = 5
   [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = 1Cu1imc2sn-e2-0-0-0-e3-0-0-0-asm-110calib
    append_date = true
  [../]
  [./perf_log]
    type =  CSV
    file_base = perf_log
    append_date = true
  [../]
  [./sim_log]
    type = Console
    append_date = true
    output_file = true
    file_base = console_sim_log
  [../]
  [./checkpoint]
    type = Checkpoint
    num_files = 2
    suffix = cp
    use_displaced = true
  [../]
[]
