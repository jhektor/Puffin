[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
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

[AuxVariables]
  # phase concentration  Sn in Cu
  # phase concentration  Sn in Cu6Sn5
  # phase concentration  Sn in Sn
  # order parameter Cu
  # order parameter Cu6Sn5
  # order parameter Sn
  [./c_cu]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.01
  [../]
  [./c_imc]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.455
  [../]
  [./c_sn]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.99
  [../]
  [./eta_cu]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta_imc]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta_sn]
    # initial_condition = 0.0
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Variables]
  # concentration Sn
  # chemical potential
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./eta1]
    # Cu
    # function = 'if(y<=10,1,0)'
    variable = eta_cu
    type = FunctionIC
    function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r<=10,1,0)'
  [../]
  [./eta2]
    # Cu6Sn5
    # function = 'if(y>10&y<=18,1,0)'
    variable = eta_imc
    type = FunctionIC
    function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>10&r<=18,1,0)'
  [../]
  [./eta3]
    # Sn
    # function = 'if(y>18,1,0)'
    variable = eta_sn
    type = FunctionIC
    function = 'r:=sqrt((x-20)^2+(y-20)^2);if(r>18,1,0)'
  [../]
  [./c]
    # Concentration of Sn
    # args = 'eta_cu eta_imc eta_sn'
    # function = '0.2*if(y<=10,1,0)+0.5*if(y>10&y<=18,1,0)+0.8*if(y>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    # function = '0.2*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.5*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.8*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    variable = c
    type = FunctionIC
    function = '0.01*if(sqrt((x-20)^2+(y-20)^2)<=10,1,0)+0.455*if(sqrt((x-20)^2+(y-20)^2)>10&sqrt((x-20)^2+(y-20)^2)<=18,1,0)+0.95*if(sqrt((x-20)^2+(y-20)^2)>18,1,0)' # TODO: Make nicer, should be possible to use values of the other variables.
  [../]
[]

[Materials]
  # Free energy
  # SwitchingFunction
  # constant properties
  [./fch_cu]
    # Chemical energy Cu phase
    type = DerivativeParsedMaterial
    f_name = fch_cu
    args = c_cu
    function = 20*(c_cu-0.01)^2
  [../]
  [./fch_imc]
    # Chemical energy Cu phase
    type = DerivativeParsedMaterial
    f_name = fch_imc
    args = c_imc
    function = 20*(c_imc-0.455)^2
  [../]
  [./fch_sn]
    # Chemical energy Sn phase
    type = DerivativeParsedMaterial
    f_name = fch_sn
    args = c_sn
    function = 20*(c_sn-0.99)^2
  [../]
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
  [./constants]
    type = GenericConstantMaterial
    prop_names = M
    prop_values = 2.7
  [../]
[]

[Kernels]
  # Kernels for split Cahn-Hilliard equation without composition gradent term(?)
  # Cahn-Hilliard Equation
  # 
  [./CHBulk]
    # Gives the residual for the concentration, dF/dc-mu
    type = KKSSplitCHCRes
    variable = c
    ca = c_imc
    cb = c_sn
    fa_name = fch_imc
    fb_name = fch_sn
    args_a = c_cu
    w = w
    h_name = h_imc
  [../]
  [./dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./ckernel]
    # Gives residual for chemical potential dcdt+M\grad(mu)
    type = SplitCHWRes
    mob_name = M
    variable = w
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       ilu            nonzero'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-11
  num_steps = 20
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

