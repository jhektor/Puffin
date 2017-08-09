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

[Variables]
  # concentration Sn
  # order parameter 1 Cu
  # order parameter 2 Cu6Sn5
  # order parameter 3 Sn
  # phase concentration 1 Sn in Cu
  # phase concentration 2 Sn in Cu6Sn5
  # phase concentration 3 Sn in Sn
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta3]
    # initial_condition = 0.0
    order = FIRST
    family = LAGRANGE
  [../]
  [./c1]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.01
  [../]
  [./c2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.455
  [../]
  [./c3]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.99
  [../]
[]

[ICs]
  [./eta1]
    # Cu
    # function = 'if(y<=10,1,0)'
    variable = eta1
    type = FunctionIC
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r<=10,1,0)'
  [../]
  [./eta2]
    # Cu6Sn5
    # function = 'if(y>10&y<=18,1,0)'
    variable = eta2
    type = FunctionIC
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r>10&r<=18,1,0)'
  [../]
  [./eta3]
    # Sn
    # function = 'if(y>18,1,0)'
    variable = eta3
    type = FunctionIC
    function = 'r:=sqrt((x-40)^2+(y-40)^2);if(r>18,1,0)'
  [../]
  [./c]
    # Concentration of Sn
    # function = '0.2*if(y<=10,1,0)+0.5*if(y>10&y<=18,1,0)+0.8*if(y>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    # function = '0.2*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.5*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.8*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' #TODO: Make nicer, should be possible to use values of the other variables.
    variable = c
    args = 'eta1 eta2 eta3'
    type = FunctionIC
    function = '0.01*if(sqrt((x-40)^2+(y-40)^2)<=10,1,0)+0.455*if(sqrt((x-40)^2+(y-40)^2)>10&sqrt((x-40)^2+(y-40)^2)<=18,1,0)+0.99*if(sqrt((x-40)^2+(y-40)^2)>18,1,0)' # TODO: Make nicer, should be possible to use values of the other variables.
  [../]
[]

