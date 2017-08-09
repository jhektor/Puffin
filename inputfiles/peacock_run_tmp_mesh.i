# 
# KKS toy problem in the non-split form
# 
# 
# This still needs finite difference preconditioning as the
# handcoded jacobians are not complete. Check out the split
# solve, which works with SMP preconditioning.
# 

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
  nz = 0
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  # order parameter
  # hydrogen concentration
  # hydrogen phase concentration (matrix)
  # hydrogen phase concentration (delta phase)
  [./eta]
    order = THIRD
    family = HERMITE
  [../]
  [./c]
    order = THIRD
    family = HERMITE
  [../]
  [./cm]
    order = THIRD
    family = HERMITE
    initial_condition = 0.0
  [../]
  [./cd]
    order = THIRD
    family = HERMITE
    initial_condition = 0.0
  [../]
[]

[Kernels]
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
    type = KKSCHBulk
    variable = c
    ca = cm
    cb = cd
    fa_name = fm
    fb_name = fd
    mob_name = 0.7
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = c
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
    kappa_name = 0.4
  [../]
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = -pc_factor_shift_type
  petsc_options_value = nonzero
  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-4
  num_steps = 1
  dt = 0.01
  dtmin = 0.01
[]

[Outputs]
  file_base = kks_example
  [./oversampling]
    type = Exodus
    refinements = 3
  [../]
[]

