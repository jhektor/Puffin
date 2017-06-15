[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 80
    ny = 80
    nz = 0
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 40
    zmin = 0
    zmax = 0
    elem_type = QUAD4
[]

[BCs]
    [./Periodic]
        [./all]
            auto_direction = 'x y'
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
        initial_condition = 0.15
    [../]

    # phase concentration  Sn in Cu6Sn5
    [./c_imc]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.55
    [../]

    # phase concentration  Sn in Sn
    [./c_sn]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.9
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
        #function = 'if(y<=10,1,0)'
        function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r<=8,1,0)'
    [../]
    [./eta2] #Cu6Sn5
        variable = eta_imc
        type = FunctionIC
        #function = 'if(y>10&y<=18,1,0)'
        function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>8&r<=16,1,0)'
    [../]
    [./eta3] #Sn
        variable = eta_sn
        type = FunctionIC
        #function = 'if(y>18,1,0)'
        function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>16,1,0)'
    [../]

    [./c] #Concentration of Sn
        variable = c
        #args = 'eta_cu eta_imc eta_sn'
        type = FunctionIC
        #function = '0.2*if(y<=10,1,0)+0.5*if(y>10&y<=18,1,0)+0.8*if(y>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
        #function = '0.2*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.5*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.8*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
        function = '0.15*if(sqrt((x-20)^2+(y-20)^2)<=8,1,0)+0.55*if(sqrt((x-20)^2+(y-20)^2)>8&sqrt((x-20)^2+(y-20)^2)<=16,1,0)+0.9*if(sqrt((x-20)^2+(y-20)^2)>16,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    [../]
[]

[Materials]
    #Free energy
    [./fch_cu] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_cu
        args = 'c_cu'
        function = '20*(c_cu-0.2)^2'
    [../]
    [./fch_imc] #Chemical energy Cu phase
        type = DerivativeParsedMaterial
        f_name = fch_imc
        args = 'c_imc'
        function = '10*(c_imc-0.5)^2'
    [../]
    [./fch_sn] #Chemical energy Sn phase
        type = DerivativeParsedMaterial
        f_name = fch_sn
        args = 'c_sn'
        function = '20*(c_sn-0.95)^2'
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
    [./g_cu]
      type = BarrierFunctionMaterial
      g_order = SIMPLE
      eta=eta_sn
      well_only = True
      function_name = g_sn
    [../]

    ## constant properties
    #[./constants]
    #    type = GenericConstantMaterial
    #    prop_names  = 'M'
    #    prop_values = '2.7'
    #[../]
    [./CHMobility]
        type = ParsedMaterial
        f_name = M
        args = 'eta_cu eta_imc eta_sn'
        constant_names = 'M_cu M_imc M_sn'
        constant_expressions = '15 15 15'
        material_property_names = 'h_cu h_imc h_sn'
        function = 'M_cu*h_cu+M_imc*h_imc+M_sn*h_sn'

    [../]
    # constant properties
    [./constants]
        type = GenericConstantMaterial
        prop_names  = 'L kappa gamma mu tgrad_corr_mult'
        prop_values = '1 0.5 0.5 1 0'
    [../]
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
        fa_name  = fch_imc
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
    [./ACInterface1] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_cu
      kappa_name = kappa
      mob_name = L
    [../]
    [./ACdfintdeta_cu]
      type = ACGrGrMulti
      variable = eta_cu
      v = 'eta_imc eta_sn'
      gamma_names = 'gamma gamma'
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
    [./ACInterface1] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_imc
      kappa_name = kappa
      mob_name = L
    [../]
    [./ACdfintdeta_imc]
      type = ACGrGrMulti
      variable = eta_imc
      v = 'eta_cu eta_sn'
      gamma_names = 'gamma gamma'
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
    [./ACInterface1] # L*kappa*grad\eta_i
      type = ACInterface
      variable = eta_sn
      kappa_name = kappa
      mob_name = L
    [../]
    [./ACdfintdeta_sn]
      type = ACGrGrMulti
      variable = eta_sn
      v = 'eta_cu eta_imc'
      gamma_names = 'gamma gamma'
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

  num_steps = 40
  dt = 0.5
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
[]
