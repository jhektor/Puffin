//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CRYSTALPLASTICITYSLIPRESISTANCEHWR_H
#define CRYSTALPLASTICITYSLIPRESISTANCEHWR_H

#include "CrystalPlasticitySlipResistance.h"

class CrystalPlasticitySlipResistanceHWR;

template <>
InputParameters validParams<CrystalPlasticitySlipResistanceHWR>();

/**
 * Phenomenological constitutive model slip resistance userobject class.
 */
class CrystalPlasticitySlipResistanceHWR : public CrystalPlasticitySlipResistance
{
public:
  CrystalPlasticitySlipResistanceHWR(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

  /// class for switching between different crystal lattice types
  static MooseEnum crystalLatticeTypeOptions();
protected:

  const MaterialProperty<std::vector<Real>> & _mat_prop_state_var;

  /// the variable to switch crystal lattice type
  MooseEnum _crystal_lattice_type;

  const Real _Q; // hardening parameter
  const Real _G0; // Lattice friction
  const Real _q; // self or latent hardening

};

#endif // CRYSTALPLASTICITYSLIPRESISTANCEHWR_H
