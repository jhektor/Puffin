/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPPLASTICENERGYMATERIAL_H
#define CPPLASTICENERGYMATERIAL_H

#include "DerivativeFunctionMaterialBase.h"

// Forward Declaration
class CPPlasticEnergyMaterial;

template <>
InputParameters validParams<CPPlasticEnergyMaterial>();

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class CPPlasticEnergyMaterial : public DerivativeFunctionMaterialBase
{
public:
  CPPlasticEnergyMaterial(const InputParameters & parameters);

  // virtual void initialSetup() override;

protected:
  virtual Real computeF() override;
  virtual Real computeDF(unsigned int i_var) override;
  virtual Real computeD2F(unsigned int i_var, unsigned int j_var) override;

  unsigned int _variable_size;
  const MaterialProperty<std::vector<Real>> & _mat_prop_state_var;
  const Real _Q;
  const Real _q;
  // const Real _s0;

  std::vector<unsigned int> _groups;
  std::vector<Real> _group_values;


};

#endif // CPPLASTICENERGYMATERIAL_H
