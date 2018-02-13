/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICENERGYMATERIALGREENPK2_H
#define ELASTICENERGYMATERIALGREENPK2_H

#include "DerivativeFunctionMaterialBase.h"

// Forward Declaration
class ElasticEnergyMaterialGreenPK2;
class RankTwoTensor;
class RankFourTensor;

template <>
InputParameters validParams<DerivativeFunctionMaterialBase>();

/**
 * Material class to compute the elastic free energy based on Green-Lagrange strain and PK2 stress and its derivatives
 */
class ElasticEnergyMaterialGreenPK2 : public DerivativeFunctionMaterialBase
{
public:
  ElasticEnergyMaterialGreenPK2(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeF() override;
  virtual Real computeDF(unsigned int i_var) override;
  virtual Real computeD2F(unsigned int i_var, unsigned int j_var) override;


  std::string _base_name;

  /// Stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;
  // std::vector<const MaterialProperty<RankTwoTensor> *> _dstress;
  // std::vector<std::vector<const MaterialProperty<RankTwoTensor> *> > _d2stress;

  ///@{ Elasticity tensor derivatives
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  std::vector<const MaterialProperty<RankFourTensor> *> _delasticity_tensor;
  std::vector<std::vector<const MaterialProperty<RankFourTensor> *>> _d2elasticity_tensor;
  ///@}

  ///@{ Strain and derivatives
  const MaterialProperty<RankTwoTensor> & _strain;
  std::vector<const MaterialProperty<RankTwoTensor> *> _dstrain;
  std::vector<std::vector<const MaterialProperty<RankTwoTensor> *>> _d2strain;
  ///@}
};

#endif // ELASTICENERGYMATERIALGREENPK2_H
