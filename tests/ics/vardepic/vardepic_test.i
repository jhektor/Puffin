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

[Variables]
    [./c]
        order = FIRST
        family = LAGRANGE
    [../]
    [./c1]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.1
    [../]
    [./c2]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.2
    [../]
    [./c3]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.3
    [../]
    [./c4]
        order = FIRST
        family = LAGRANGE
        initial_condition = 0.4
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
        order = FIRST
        family = LAGRANGE
    [../]
    [./eta4]
        order = FIRST
        family = LAGRANGE
    [../]
[]
[ICs]
    [./eta1]
        type = MultiBoundingBoxIC
        variable = eta1
        corners = '0. 0. 0.'
        opposite_corners = '2000. 600. 0.'
        inside = 1.
        outside = 0.
    [../]
    [./eta2]
        type = MultiBoundingBoxIC
        variable = eta2
        corners = '0. 650. 0.   1550. 650. 0.'
        opposite_corners = '500. 1000. 0   2000. 1000. 0.'
        inside = 1.
        outside = 0.
    [../]
    [./eta3]
        type = MultiBoundingBoxIC
        variable = eta3
        corners = '550. 650. 0.'
        opposite_corners = '1500. 1000. 0.'
        inside = 1.
        outside = 0.
    [../]
    [./eta4]
        type = MultiBoundingBoxIC
        variable = eta4
        corners = '0. 1050. 0.'
        opposite_corners = '2000. 3000. 0.'
        inside = 1.
        outside = 0.
    [../]
    [./c]
        type = VarDepIC
        variable = c
        etas = 'eta1 eta2 eta3 eta4'
        cis = 'c1 c2 c3 c4'
    [../]
[]
[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]
[Outputs]
  file_base = out
  exodus = true
[]
