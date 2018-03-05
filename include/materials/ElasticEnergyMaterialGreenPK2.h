/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELASTICENERGYMATERIALGREENPK2_H
#define ELASTICENERGYMATERIALGREENPK2_H

#include "DerivativeFunctionMaterialBase.h"
#include "DerivativeFunctionMaterialBase.h"
#include "RankTwoTensor.h"

// Forward Declaration
class ElasticEnergyMaterialGreenPK2;
// class RankTwoTensor;
class RankFourTensor;

template <>
InputParameters validParams<DerivativeFunctionMaterialBase>();

/**
 * Material class to compute the elastic free energy based on Green-Lagrange strain and PK2 stress
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


  std::string _base_name, _eigenstrain_name;
  bool _plasticity, _eigenstrain;


  // RankTwoTensor _firr;
  const MaterialProperty<RankTwoTensor> & _f;
  const MaterialProperty<RankTwoTensor> & _firr;

  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  const MaterialProperty<RankTwoTensor> & _stress;



};

#endif // ELASTICENERGYMATERIALGREENPK2_H
