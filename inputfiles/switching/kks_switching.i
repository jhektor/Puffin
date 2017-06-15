#
# This test is for the 3-phase KKS model
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  nz = 0
  xmin = 0
  xmax = 80
  ymin = 0
  ymax = 80
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

  # order parameter 1 Cu
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter 2 Cu6Sn5
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter 3 Sn
  [./eta3]
    order = FIRST
    family = LAGRANGE
    #initial_condition = 0.0
  [../]

  # phase concentration 1 Sn in Cu
  [./c1]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.01
  [../]

  # phase concentration 2 Sn in Cu6Sn5
  [./c2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.455
  [../]

  # phase concentration 3 Sn in Sn
  [./c3]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.99
  [../]
[]

[ICs]
  [./eta1] #Cu
    variable = eta1
    type = FunctionIC
    #function = 'if(y<=10,1,0)'
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r<=10,1,0)'
  [../]
  [./eta2] #Cu6Sn5
    variable = eta2
    type = FunctionIC
    #function = 'if(y>10&y<=18,1,0)'
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r>10&r<=18,1,0)'

  [../]
  [./eta3] #Sn
    variable = eta3
    type = FunctionIC
    #function = 'if(y>18,1,0)'
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r>18,1,0)'

  [../]
  [./c] #Concentration of Sn
    variable = c
    args = 'eta1 eta2 eta3'
    type = FunctionIC
    #function = '0.2*if(y<=10,1,0)+0.5*if(y>10&y<=18,1,0)+0.8*if(y>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    #function = '0.2*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.5*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.8*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    function = '0.01*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.455*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.99*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
  [../]
[]


[Materials]
  # simple toy free energies
  [./f1]
    type = DerivativeParsedMaterial
    f_name = F1
    args = 'c1'
    #function = '20*(c1-0.2)^2'
    function = '20*(c1-0.01)^2'
  [../]
  [./f2]
    type = DerivativeParsedMaterial
    f_name = F2
    args = 'c2'
    #function = '20*(c2-0.5)^2'
    function = '20*(c2-0.455)^2'
  [../]
  [./f3]
    type = DerivativeParsedMaterial
    f_name = F3
    args = 'c3'
    #function = '20*(c3-0.8)^2'
    function = '20*(c3-0.99)^2'
  [../]

  # Switching functions for each phase
  # h1(eta1, eta2, eta3)
  [./h1]
    type = SwitchingFunctionMultiPhaseMaterial
    all_etas = 'eta1 eta2 eta3'
    phase_etas = eta1
    h_name = h1
  [../]
  # h2(eta1, eta2, eta3)
  [./h2]
    type = SwitchingFunctionMultiPhaseMaterial
    all_etas = 'eta1 eta2 eta3'
    phase_etas = eta2
    h_name = h2
  [../]
  # h3(eta1, eta2, eta3)
  [./h3]
    type = SwitchingFunctionMultiPhaseMaterial
    all_etas = 'eta1 eta2 eta3'
    phase_etas = eta3
    h_name = h3
  [../]



[Kernels]
  # Kernels for Allen-Cahn equation for eta1
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]



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

  num_steps = 200
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
