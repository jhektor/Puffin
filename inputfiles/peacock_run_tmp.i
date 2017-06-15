# 
# KKS toy problem in the split form
# 
# 
# Precondition using handcoded off-diagonal terms
# 

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 15
  ny = 15
  nz = 0
  xmin = -2.5
  xmax = 2.5
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[AuxVariables]
  [./Fglobal]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Variables]
  # order parameter
  # hydrogen concentration
  # chemical potential
  # hydrogen phase concentration (matrix)
  # hydrogen phase concentration (delta phase)
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
  [./cm]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./cd]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
[]

[ICs]
  [./eta]
    variable = eta
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 1.5
    invalue = 0.2
    outvalue = 0.1
    int_width = 0.75
  [../]
  [./c]
    variable = c
    type = SmoothCircleIC
    x1 = 0.0
    y1 = 0.0
    radius = 1.5
    invalue = 0.6
    outvalue = 0.4
    int_width = 0.75
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = 'eta w c cm cd'
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  # Free energy of the matrix
  # Free energy of the delta phase
  # h(eta)
  # g(eta)
  # constant properties
  [./fm]
    type = DerivativeParsedMaterial
    f_name = fm
    args = cm
    function = (0.1-cm)^2
  [../]
  [./fd]
    type = DerivativeParsedMaterial
    f_name = fd
    args = cd
    function = (0.9-cd)^2
  [../]
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
  [../]
  [./g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  [../]
  [./constants]
    type = GenericConstantMaterial
    prop_names = 'M   L   kappa'
    prop_values = '0.7 0.7 0.4  '
  [../]
[]

[Kernels]
  # full transient
  # enforce c = (1-h(eta))*cm + h(eta)*cd
  # enforce pointwise equality of chemical potentials
  # 
  # Cahn-Hilliard Equation
  # 
  # 
  # Allen-Cahn Equation
  # 
  [./PhaseConc]
    type = KKSPhaseConcentration
    ca = cm
    variable = cd
    c = c
    eta = eta
  [../]
  [./ChemPotVacancies]
    type = KKSPhaseChemicalPotential
    variable = cm
    cb = cd
    fa_name = fm
    fb_name = fd
  [../]
  [./CHBulk]
    type = KKSSplitCHCRes
    variable = c
    ca = cm
    cb = cd
    fa_name = fm
    fb_name = fd
    w = w
  [../]
  [./dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./ckernel]
    type = SplitCHWRes
    mob_name = M
    variable = w
  [../]
  [./ACBulkF]
    type = KKSACBulkF
    variable = eta
    fa_name = fm
    fb_name = fd
    args = 'cm cd'
    w = 0.4
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = eta
    ca = cm
    cb = cd
    fa_name = fm
    fb_name = fd
  [../]
  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  [../]
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
[]

[AuxKernels]
  [./GlobalFreeEnergy]
    variable = Fglobal
    type = KKSGlobalFreeEnergy
    fa_name = fm
    fb_name = fd
    w = 0.4
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = -pc_factor_shift_type
  petsc_options_value = nonzero
  l_max_its = 100
  nl_max_its = 100
  num_steps = 3
  dt = 0.1
[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Outputs]
  file_base = kks_example_split
  exodus = true
[]

