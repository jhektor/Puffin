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
    [./c]
        order = FIRST
        family = LAGRANGE
    [../]
    # chemical potential
    [./w]
        order = FIRST
        family = LAGRANGE
    [../]

    # phase concentration  Sn in Cu
    [./c_cu]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.05
    [../]

    # phase concentration  Sn in Cu6Sn5
    [./c_imc]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.455
    [../]

    # phase concentration  Sn in Sn
    [./c_sn]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.95 #0.95#0.999
    [../]

    # order parameter Cu
    [./eta_cu]
        order = FIRST
        family = LAGRANGE
    [../]
    # order parameter Cu6Sn5
    [./eta_imc]
        order = FIRST
        family = LAGRANGE
    [../]

    # order parameter Sn
    [./eta_sn]
        order = FIRST
        family = LAGRANGE
    [../]

[]
[ICs]
    [./eta1] #Cu
        variable = eta_cu
        type = FunctionIC
        function = 'if(x<=48000,1,0)'
    [../]
    [./eta2] #Cu6Sn5
        variable = eta_imc
        type = FunctionIC
        function = 'if(x>48000&x<=52000,1,0)'
    [../]
    [./eta3] #Sn
        variable = eta_sn
        type = FunctionIC
        function = 'if(x>52000,1,0)'
    [../]

    [./c] #Concentration of Sn
        variable = c
        type = FunctionIC
        function = 'if(x<=48000,0.05,0)+if(x>48000&x<=52000,0.455,0)+if(x>52000,0.95,0)'
    [../]
[]

[Materials]
    #scalings
    [./scale]
        type = GenericConstantMaterial
        prop_names = 'length_scale energy_scale time_scale'
        prop_values = '1e9 6.24150943e18 1.' #m to nm J to eV s to h
    [../]
    #Constants
    [./energy_constants]
        type = GenericConstantMaterial
        prop_names = 'A_cu A_imc A_sn c_cu_imc c_imc_cu c_imc_sn c_sn_imc c_cu_sn c_sn_cu'
        prop_values = '1e8 1e9 1e9 0.0504 0.4524 0.4760 0.999 0.1938 0.9898' #J/m^3  -
        #prop_values = '1e9 1e9 1e9 0.0504 0.4524 0.4760 0.999 0.1938 0.9898' #J/m^3  -
    [../]
    [./diffusion_constants]
      type = GenericConstantMaterial
      prop_names = 'D_cu D_imc D_sn'
      prop_values = '1e-25 1e-16 1e-14' # m^2/s
      #prop_values = '1e-25 1e-16 1e-13' # m^2/s
      #prop_values = '1e-17 1e-16 1e-16' # m^2/s
    [../]
    [./model_constants]
      type = GenericConstantMaterial
      prop_names = 'sigma delta gamma tgrad_corr_mult'
      prop_values = '0.5 0.666e-6 1.5 0' #J/m^2 m - ?
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
    #Free energy
    [./fch_cu] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_cu
        args = 'c_cu'
        material_property_names = 'A_cu length_scale energy_scale'
        function = '(energy_scale/length_scale^3)*(0.5*A_cu*(c_cu-0.076)^2)' #eV/nm^3
        derivative_order = 2
    [../]
    [./fch_imc] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_imc
        args = 'c_imc'
        material_property_names = 'A_imc length_scale energy_scale'
        function = '(energy_scale/length_scale^3)*(0.5*A_imc*(c_imc-0.455)^2-1e6)' #eV/nm^3
        derivative_order = 2
    [../]
    [./fch_sn] #Chemical energy Sn phase
        type = DerivativeParsedMaterial
        f_name = fch_sn
        args = 'c_sn'
        material_property_names = 'A_sn length_scale energy_scale'
        function = '(energy_scale/length_scale^3)*(0.5*A_sn*(c_sn-0.978)^2+1e7)' #eV/nm^3
        derivative_order = 2
    [../]

    #SwitchingFunction
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

    #Double well, not used MAYBE USE TO KEEP THE ORDER PARAMETERS IN [0:1]
    [./g_cu]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_cu
      well_only = True
      function_name = g_cu
    [../]
    #Double well, not used
    [./g_imc]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_imc
      well_only = True
      function_name = g_imc

    [../]
    #Double well, not used
    [./g_sn]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_sn
      well_only = True
      function_name = g_sn
    [../]
    [./CHMobility]
        type = DerivativeParsedMaterial
        f_name = M
        args = 'eta_cu eta_imc eta_sn'
        material_property_names = 'h_cu(eta_cu,eta_imc,eta_sn) h_imc(eta_cu,eta_imc,eta_sn) h_sn(eta_cu,eta_imc,eta_sn) D_cu D_imc D_sn A_cu A_imc A_sn length_scale energy_scale time_scale'
        function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_cu/A_cu+h_imc*D_imc/A_imc+h_sn*D_sn/A_sn)' #nm^5/eVs
        #function = '(length_scale^5/(energy_scale*time_scale))*(h_cu*D_sn/A_sn+h_imc*D_sn/A_sn+h_sn*D_sn/A_sn)' #nm^5/eVs
        derivative_order = 2
    [../]

    [./ACMobility]
        type = DerivativeParsedMaterial
        f_name = L
        args = 'eta_cu eta_imc eta_sn'
        material_property_names = 'L_cu_imc L_imc_sn L_cu_sn' # h_cu(eta_cu,eta_imc,eta_sn) h_imc(eta_cu,eta_imc,eta_sn) h_sn(eta_cu,eta_imc,eta_sn)'
        #function ='(L_cu_imc*eta_cu^2*eta_imc^2+L_imc_sn*eta_imc^2*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2)/(eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2)'

        #function ='0.5*(L_cu_imc+L_imc_sn)' #L_imc_sn'

        # Added epsilon to prevent division by 0 (Larry Aagesen) 0 i bulk funkar inte
        function ='pf:=1e5;eps:=0.01;if(eta_cu>0.99|eta_imc>0.99|eta_sn>0.99,0,(L_cu_imc*(pf*eta_cu^2+eps)*(pf*eta_imc^2+eps)+L_imc_sn*(pf*eta_imc^2+eps)*(pf*eta_sn^2+eps)+L_cu_sn*(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps))/((pf*eta_cu^2+eps)*(pf*eta_imc^2+eps)+(pf*eta_imc^2+eps)*(pf*eta_sn^2+eps)+(pf*eta_cu^2+eps)*(pf*eta_sn^2+eps)))'
        #function ='pf:=1e2;eps:=0.01;(L_cu_imc*(pf*eta_cu+eps)^2*(pf*eta_imc+eps)^2+L_imc_sn*(pf*eta_imc+eps)^2*(pf*eta_sn+eps)^2+L_cu_sn*(pf*eta_cu+eps)^2*(pf*eta_sn+eps)^2)/((pf*eta_cu+eps)^2*(pf*eta_imc+eps)^2+(pf*eta_imc+eps)^2*(pf*eta_sn+eps)^2+(pf*eta_cu+eps)^2*(pf*eta_sn+eps)^2)'

        # Conditional function (Daniel Schwen)
        #function ='numer:=L_cu_imc*eta_cu^2*eta_imc^2+L_imc_sn*eta_imc^2*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2;denom:=eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2;if(denom!=0,numer/denom,0.5*(L_cu_imc+L_imc_sn))'
        #function ='numer:=L_cu_imc*eta_cu^2*eta_imc^2+L_imc_sn*eta_imc^2*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2;denom:=eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2;if(denom!=0,100*numer/denom,0.5*(L_cu_imc+L_imc_sn))'
        #function ='numer:=L_cu_imc*eta_cu^2*eta_imc^2+L_imc_sn*eta_imc^2*eta_sn^2+L_cu_sn*eta_cu^2*eta_sn^2;denom:=eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2;if(denom!=0,numer/denom,0)'

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
        ca       = c_imc
        cb       = c_sn
        fa_name  = fch_imc #only fa is used
        fb_name  = fch_sn
        w        = w
        h_name   = h_imc
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
        args = 'eta_cu eta_imc eta_sn'
    [../]

    #KKS conditions
    # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n
    [./chempot_cu_imc]
      type = KKSPhaseChemicalPotential
      variable = c_cu
      cb       = c_imc
      fa_name  = fch_cu
      fb_name  = fch_imc
    [../]
    [./chempot_sn_cu]
      type = KKSPhaseChemicalPotential
      variable = c_imc
      cb       = c_sn
      fa_name  = fch_imc
      fb_name  = fch_sn
    [../]
    [./phaseconcentration] # enforce c = sum h_i*c_i
      type = KKSMultiPhaseConcentration
      variable = c_sn
      cj = 'c_cu c_imc c_sn'
      hj_names = 'h_cu h_imc h_sn'
      etas = 'eta_cu eta_imc eta_sn'
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
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      gi_name   = g_cu
      eta_i     = eta_cu
      wi        = 5
      mob_name = L
      args      = 'c_cu c_imc c_sn eta_imc eta_sn'
    [../]
    [./ACBulkC_cu] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_cu
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      cj_names  = 'c_cu c_imc c_sn'
      eta_i     = eta_cu
      mob_name = L
      args      = 'eta_imc eta_sn'
    [../]
    [./ACInterface_cu] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_cu
      kappa_name = kappa
      mob_name = L
      args      = 'eta_imc eta_sn'
      variable_L = true
    [../]
    [./ACdfintdeta_cu] #L*m*(eta_i^3-eta_i+2*beta*eta_i*sum_j eta_j^2)
      type = ACGrGrMulti
      variable = eta_cu
      v = 'eta_imc eta_sn'
      gamma_names = 'gamma gamma'
      mob_name = L
      args = 'eta_imc eta_sn'
    [../]

    #Kernels for Allen-Cahn equation for Cu6Sn5
    [./detadt_imc]
      type = TimeDerivative
      variable = eta_imc
    [../]
    [./ACBulkF_imc] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_imc
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      gi_name   = g_imc
      eta_i     = eta_imc
      wi        = 5
      mob_name = L
      args      = 'c_cu c_imc c_sn eta_cu eta_sn'
    [../]
    [./ACBulkC_imc] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_imc
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      cj_names  = 'c_cu c_imc c_sn'
      eta_i     = eta_imc
      mob_name = L
      args      = 'eta_cu eta_sn'
    [../]
    [./ACInterface_imc] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_imc
      kappa_name = kappa
      mob_name = L
      args      = 'eta_cu eta_sn'
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

    #Kernels for Allen-Cahn equation for Sn
    [./detadt_sn]
      type = TimeDerivative
      variable = eta_sn
    [../]
    [./ACBulkF_sn] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_sn
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      gi_name   = g_sn
      eta_i     = eta_sn
      wi        = 5
      mob_name = L
      args      = 'c_cu c_imc c_sn eta_imc eta_cu'
    [../]
    [./ACBulkC_sn] # -L\sum_j dh_j/deta_i*mu_jc_j
      type = KKSMultiACBulkC
      variable  = eta_sn
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      cj_names  = 'c_cu c_imc c_sn'
      eta_i     = eta_sn
      mob_name = L
      args      = 'eta_cu eta_imc'
    [../]
    [./ACInterface_sn] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_sn
      kappa_name = kappa
      mob_name = L
      args      = 'eta_cu eta_imc'
      variable_L = true
    [../]
    [./ACdfintdeta_sn]
      type = ACGrGrMulti
      variable = eta_sn
      v = 'eta_cu eta_imc'
      gamma_names = 'gamma gamma'
      mob_name = L
      args= 'eta_cu eta_imc'
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
        hj_names = 'h_cu h_imc h_sn'
        Fj_names = 'fch_cu fch_imc fch_sn'
        gj_names = 'g_cu g_imc g_sn'
        additional_free_energy = f_int
        interfacial_vars = 'eta_cu eta_imc eta_sn'
        kappa_names = 'kappa kappa kappa'
        w = 5
        execute_on = 'initial timestep_end'
    [../]
    [./f_int]
        type = ParsedAux
        variable = f_int
        args = 'eta_cu eta_imc eta_sn'
        constant_names = 'sigma delta gamma length_scale energy_scale'
        constant_expressions = '0.5 0.667e-6 1.5 1e9 6.24150943e18'
        function ='mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(0.25*eta_cu^4-0.5*eta_cu^2+0.25*eta_imc^4-0.5*eta_imc^2+0.25*eta_sn^4-0.5*eta_sn^2+gamma*(eta_cu^2*(eta_imc^2+eta_sn^2)+eta_imc^2*eta_sn^2)+0.25)'
        execute_on = 'initial timestep_end'
    [../]
    [./s]
      type = ParsedAux
      variable = s
      args = 'eta_cu eta_imc eta_sn'
      #function = 'eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2'
      function = 'eta_cu+eta_imc+eta_sn'
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
    [./imc_area_h]
      type = ElementIntegralMaterialProperty
      mat_prop = h_imc
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
  show_var_residual_norms = false
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
  nl_rel_tol = 1.0e-9 #1.0e-10
  nl_abs_tol = 1.0e-10#1.0e-11

  #num_steps = 2000
  end_time = 64e6
  #very simple adaptive time stepper
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
  file_base = keep/moelans2011fig2_Lvar0bulk_Mvar_sharp_barrier5_100um
  [./exodus_out]
    type = Exodus
    interval = 1
  [../]
  #exodus = true
  csv = true
  print_perf_log = true
  interval = 1 #5
[]
