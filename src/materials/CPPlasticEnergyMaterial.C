/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPPlasticEnergyMaterial.h"
// #include "CrystalPlasticityStateVariable.h"

template <>
InputParameters
validParams<CPPlasticEnergyMaterial>()
{
  InputParameters params = validParams<DerivativeFunctionMaterialBase>();
  params.addClassDescription("Free energy material for the plastic energy contributions. Based on Mellbin 2017");
  params.addRequiredParam<std::string>("uo_state_var_name",
                               "Name of state variable property: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addRequiredCoupledVar("args", "Arguments of F() - use vector coupling");
  params.addParam<Real>("Q", 1, "Hardening parameter SHOULD BE ALREADY USED SOMEWHERE...");
  params.addParam<Real>("q", 1.4, "Cross-hardening parameter");
  params.addParam<Real>("s0", 0, "Initial value of slip resistance"); //Should make a vector
  params.addRequiredParam<unsigned int>("variable_size", "Number of slip systems.");

  return params;
}

CPPlasticEnergyMaterial::CPPlasticEnergyMaterial(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
  _variable_size(getParam<unsigned int>("variable_size")),
  _mat_prop_state_var(
    getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
  _Q(getParam<Real>("Q")),
  _q(getParam<Real>("q")),
  _s0(getParam<Real>("s0"))
{

}
Real
CPPlasticEnergyMaterial::computeF()
{
  Real sum, qab;
  sum = 0;
  // Loop over all slip systems
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < _variable_size; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      // Not sure why not if i==j
      // if (iplane == jplane) // Kalidindi
      if (i == j)
        qab = 1.0;
      else
        qab = _q;
      sum += qab * _mat_prop_state_var[_qp][i] * _mat_prop_state_var[_qp][j];
      // This will not give sum = 0 for purely elastic so must subtract initial value of _mat_prop_state_var
      sum -= qab * _s0 * _s0;
    }
  }

  return 0.5 * _Q * sum;
}

/// Do not know which derivatives I must calculate so just return 0
Real
CPPlasticEnergyMaterial::computeDF(unsigned int i_var)
{
  return 0;
}

Real
CPPlasticEnergyMaterial::computeD2F(unsigned int i_var, unsigned int j_var)
{
  return 0;
}
