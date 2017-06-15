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

[Outputs]
  exodus = true
[]

