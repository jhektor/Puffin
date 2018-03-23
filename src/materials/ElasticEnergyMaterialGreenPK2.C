/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElasticEnergyMaterialGreenPK2.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

template <>
InputParameters
validParams<ElasticEnergyMaterialGreenPK2>()
{
  InputParameters params = validParams<DerivativeFunctionMaterialBase>();
  params.addClassDescription("Free energy material for the elastic energy contributions. This calculates the energy as 0.5*S:E. THIS DOES NOT CALCULATE ANY DERIVATIVES.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<std::string>("eigenstrain_name", "Name of eigenstrain material");
  params.addParam<bool>("plasticity",false,"Plasticity in this phase?");
  params.addParam<bool>("eigenstrain",false,"Eigenstrain in this phase?");
  params.addRequiredCoupledVar("args", "Arguments of F() - use vector coupling");
  params.addCoupledVar("displacement_gradients",
                       "Vector of displacement gradient variables (see "
                       "Modules/PhaseField/DisplacementGradients "
                       "action)");

  return params;
}

ElasticEnergyMaterialGreenPK2::ElasticEnergyMaterialGreenPK2(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain_name(isParamValid("eigenstrain_name") ? _base_name + getParam<std::string>("eigenstrain_name") : ""),
    _plasticity(getParam<bool>("plasticity")),
    _eigenstrain(getParam<bool>("eigenstrain")),
    _f(getMaterialProperty<RankTwoTensor>(_base_name  +"deformation_gradient")),
    _firr(_plasticity ? getMaterialPropertyByName<RankTwoTensor>(_base_name  +"fp") : (_eigenstrain ? getMaterialPropertyByName<RankTwoTensor>(_eigenstrain_name) : _f)),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _stress(_plasticity ? getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2") : getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress"))
{


}

void
ElasticEnergyMaterialGreenPK2::initialSetup()
{
}

Real
ElasticEnergyMaterialGreenPK2::computeF()
{
  RankTwoTensor fe, E, S, I, fstar;
  I.zero();
  I.addIa(1.0);

  if (_plasticity)
    fe = _f[_qp] * _firr[_qp].inverse(); // Multiplicative split of deformation gradient
  else if (_eigenstrain)
  {
    fstar = _firr[_qp] + I; // Check if this makes sense
    fe = _f[_qp] * fstar.inverse();
  }
  else
    fe = _f[_qp];

  // Elastic Green-Lagrange Strain
  E = 0.5 * (fe.transpose() * fe -  I);
  // second Piola-Kirchoff stress
  if (_plasticity)
    S = _stress[_qp]; // this is pk2 from Crystal Plasticity, with rotated elasticity_tensor
  else
    S = _elasticity_tensor[_qp] * E;

  return 0.5 * S.doubleContraction(E);
}

Real
ElasticEnergyMaterialGreenPK2::computeDF(unsigned int i_var)
{
  return 0.0;
}

Real
ElasticEnergyMaterialGreenPK2::computeD2F(unsigned int i_var, unsigned int j_var)
{
  return 0.0;
}
