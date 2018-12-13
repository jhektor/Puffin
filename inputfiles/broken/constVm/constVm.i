[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 100
    ny = 100
    nz = 0
    xmin = 0
    xmax = 3
    ymin = 0
    ymax = 3
    zmin = 0
    zmax = 0
    elem_type = QUAD4
[]

[BCs]
    [./Periodic]
        [./x]
            variable = 'c c_cu c_imc c_sn w eta_cu eta_imc eta_sn'
            auto_direction = 'x'
        [../]
    [../]
    [./nonconserved]
        type = NeumannBC
        boundary = 'bottom top'
        variable = 'c c_cu c_imc c_sn w'
    [../]
    [./conserved]
        type = NeumannBC
        variable = 'eta_cu eta_imc eta_sn'
        boundary = 'bottom top'
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
        initial_condition = 0.01
    [../]

    # phase concentration  Sn in Cu6Sn5
    [./c_imc]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.417
    [../]

    # phase concentration  Sn in Sn
    [./c_sn]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.99
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
    #initial_condition = 0.0
    [../]


[]

[ICs]
    [./eta1] #Cu
        variable = eta_cu
        type = FunctionIC
        #function = 'if(y<=0.8,1,0)'
        function = 'r:=sqrt((x-1.5)^2+(y-1.5)^2);if(r<=0.25,1,0)'
    [../]
    [./eta2] #Cu6Sn5
        variable = eta_imc
        type = FunctionIC
        #function = 'if(y>0.8&y<=1.3,1,0)'
        function = 'r:=sqrt((x-1.5)^2+(y-1.5)^2);if(r>0.25&r<=0.5,1,0)'
    [../]
    [./eta3] #Sn
        variable = eta_sn
        type = FunctionIC
        #function = 'if(y>1.3,1,0)'
        function = 'r:=sqrt((x-1.5)^2+(y-1.5)^2);if(r>0.5,1,0)'
    [../]

    [./c] #Concentration of Sn
        variable = c
        #args = 'eta_cu eta_imc eta_sn'
        type = FunctionIC
        #function = '(0.01*if(y<=0.8,1,0)+0.417*if(y>0.8&y<=1.3,1,0)+0.99*if(y>1.3,1,0))' #[-] #TODO: Make nicer, should be possible to use values of the other variables.
        function = '0.01*if(sqrt((x-1.5)^2+(y-1.5)^2)<=0.25,1,0)+0.417*if(sqrt((x-1.5)^2+(y-1.5)^2)>0.25&sqrt((x-1.5)^2+(y-1.5)^2)<=0.5,1,0)+0.99*if(sqrt((x-1.5)^2+(y-1.5)^2)>0.5,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
        #function = '0.01*if(sqrt((x-20)^2+(y-20)^2)<=8,1,0)+0.417*if(sqrt((x-20)^2+(y-20)^2)>8&sqrt((x-20)^2+(y-20)^2)<=16,1,0)+0.99*if(sqrt((x-20)^2+(y-20)^2)>16,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    [../]
[]

[Materials]
    # Constants (temp = 298 K)
    [./energy_A]
        type = GenericConstantMaterial
        prop_names = 'A_cu A_imc A_sn'
        prop_values = '1.0133e5 4e5 4.2059e6' #[J/mol]
    [../]
    [./energy_B]
        type = GenericConstantMaterial
        prop_names = 'B_cu B_imc B_sn'
        prop_values = '-2.1146e4 -6.9892e3 7.1680e3' #[J/mol]
    [../]
    [./energy_C]
        type = GenericConstantMaterial
        prop_names = 'Cc_cu Cc_imc Ccsn'
        prop_values = '-1.2842e4 -1.9185e4 -1.5265e4' #[J/mol]
    [../]
    [./energy_xhat]
        type = GenericConstantMaterial
        prop_names = 'chat_cu chat_imc chat_sn'
        prop_values = '0.10569 0.41753 0.99941' #[J/mol]
    [../]
    [./xhat_eq]
        type = GenericConstantMaterial
        prop_names = 'chat_eq_cu chat_eq_imc_cu chat_eq_imc_sn chat_eq_sn'
        prop_values = '0.10569 0.38210 0.45290 0.99941' #[J/mol]
    [../]
    [./molar_volume]
        type = GenericConstantMaterial
        prop_names = 'Vm'
        prop_values = '16.29e12' #[um^3/mol]
    [../]
    [./diffusion_coeff]
        type = GenericConstantMaterial
        prop_names = 'D_cu D_imc D_sn'
        prop_values = '2.877e-24 6.575e-7 2.452e-17' #[um^2/s]
    [../]
    [./model_parameters]
        type = GenericConstantMaterial
        prop_names  = 'sigma delta tgrad_corr_mult' #surface energy[J/um^2] interface width[um] ??
        prop_values = '0.5e-12 6.*dx 0'
    [../]
    [./derived_parameters]
        type = GenericConstantMaterial
        prop_names = 'mu alpha beta'
        prop_values = '6.*sigma/delta 3.*delta*sigma/4. 1.5'
    [../]

    #Free energy
    [./fch_cu] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_cu
        args = 'c_cu'
        material_property_names = 'A_cu B_cu Cc_cu chat_cu Vm'
        function = '(0.5*A_cu*(c_cu-chat_cu)^2+B_cu*(c_cu-chat_cu)+Cc_cu)/Vm' #[J/um^3]
    [../]
    [./fch_imc] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_imc
        args = 'c_imc'
        material_property_names = 'A_imc B_imc Cc_imc chat_imc Vm'
        function = '(0.5*A_imc*(c_imc-chat_imc)^2+B_imc*(c_imc-chat_imc)+Cc_imc)/Vm'
    [../]
    [./fch_sn] #Chemical energy Sn phase
        type = DerivativeParsedMaterial
        f_name = fch_sn
        args = 'c_sn'
        material_property_names = 'A_sn B_sn Cc_sn chat_sn Vm'
        function = '(0.5*A_sn*(c_sn-chat_sn)^2+B_sn*(c_sn-chat_sn)+Cc_sn)/Vm'
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

    #Double well, not used
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

    ## constant properties
    #[./constants]
    #    type = GenericConstantMaterial
    #    prop_names  = 'M L'
    #    prop_values = '1e10. 1.'
    #[../]
    [./CHMobility]
        type = DerivativeParsedMaterial
        f_name = M
        args = 'eta_cu eta_imc eta_sn'
        material_property_names = 'h_cu h_imc h_sn D_cu D_imc D_sn Vm A_cu A_imc A_sn'
        function = '(h_cu*D_cu/A_cu/Vm+h_imc*D_imc/A_imc/Vm+h_sn*D_sn/A_sn/Vm)'
        outputs = exodus
    [../]
    [./ACMobility] #not correct
        type = DerivativeParsedMaterial
        f_name = L
        args = 'eta_cu eta_imc eta_sn'
        material_property_names = 'M'
        function = '100.*M'
        outputs = exodus
    [../]
    #[./ACMobility] #not correct
    #    type = DerivativeParsedMaterial
    #    f_name = L
    #    args = 'eta_cu eta_imc eta_sn'
    #    material_property_names = 'mu Vm D_cu A_cu D_imc A_imc D_sn A_sn alpha chat_eq_cu chat_eq_imc_cu chat_eq_imc_sn chat_eq_sn'
    #    function='(eta_cu^2*eta_imc^2*mu*Vm*(D_cu/A_cu+D_imc/A_imc)/(3.*alpha*(chat_eq_cu-chat_eq_imc_cu)^2)+eta_imc^2*eta_sn^2*mu*Vm*(D_sn/A_sn+D_imc/A_imc)/(3.*alpha*(chat_eq_sn-chat_eq_imc_sn)^2)+eta_cu^2*eta_sn^2*mu*Vm*(D_cu/A_cu+D_sn/A_sn)/(3.*alpha*(chat_eq_cu-chat_eq_sn)^2))/(eta_cu^2*eta_imc^2+eta_imc^2*eta_sn^2+eta_cu^2*eta_sn^2)'
    #[../]
[]

[Kernels]
    #Kernels for split Cahn-Hilliard equation without composition gradent term(?)
    # Cahn-Hilliard Equation
    #
    [./CHBulk] # Gives the residual for the concentration, dF/dc-mu
        type = KKSSplitCHCRes
        variable = c
        ca       = c_imc
        cb       = c_sn
        fa_name  = fch_imc #only fa is used
        fb_name  = fch_sn
        #args_a = 'c_cu'
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
    [./ACBulkF_cu] # sum_j dh_j/deta_i*F_j+w*dg/deta_i, last term is not used
      type = KKSMultiACBulkF
      variable  = eta_cu
      Fj_names  = 'fch_cu fch_imc fch_sn'
      hj_names  = 'h_cu h_imc h_sn'
      gi_name   = g_cu
      eta_i     = eta_cu
      wi        = 0
      mob_name = L
      args      = 'eta_imc eta_sn c_cu c_imc c_sn'
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
      kappa_name = alpha
      mob_name = L
      args      = 'eta_imc eta_sn'
    [../]
    [./ACdfintdeta_cu]
      type = ACGrGrMulti
      variable = eta_cu
      v = 'eta_imc eta_sn'
      gamma_names = 'beta beta'
      mob_name = L
      args      = 'eta_imc eta_sn'
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
      wi        = 0
      mob_name = L
      args      = 'eta_cu eta_sn c_cu c_imc c_sn'
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
      kappa_name = alpha
      mob_name = L
      args = 'eta_cu eta_sn'
    [../]
    [./ACdfintdeta_imc]
      type = ACGrGrMulti
      variable = eta_imc
      v = 'eta_cu eta_sn'
      gamma_names = 'beta beta'
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
      wi        = 0
      mob_name = L
      args      = 'eta_cu eta_imc c_cu c_imc c_sn'
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
      kappa_name = alpha
      mob_name = L
      args = 'eta_cu eta_imc'
    [../]
    [./ACdfintdeta_sn]
      type = ACGrGrMulti
      variable = eta_sn
      v = 'eta_cu eta_imc'
      gamma_names = 'beta beta'
      mob_name = L
      args = 'eta_cu eta_imc'
    [../]

[]

[Postprocessors]
    [./imc_thickness]
        type = NodalSum
        execute_on = timestep_end
        variable = eta_imc
    [../]
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-11

  num_steps = 100
  dt = 250
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
