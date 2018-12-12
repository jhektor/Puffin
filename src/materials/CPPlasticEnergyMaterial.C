/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPPlasticEnergyMaterial.h"
// #include "CrystalPlasticityStateVariable.h"
registerMooseObject("PuffinApp", CPPlasticEnergyMaterial);

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
  // params.addParam<Real>("s0", 0, "Initial value of slip resistance"); //Should make a vector
  params.addRequiredParam<unsigned int>("variable_size", "Number of slip systems.");
  params.addParam<std::vector<unsigned int>>("groups",
                                             "To group the initial values on different "
                                             "slip systems 'format: [start end)', i.e.'0 "
                                             "4 8 11' groups 0-3, 4-7 and 8-11 ");
  params.addParam<std::vector<Real>>("group_values",
                                     "The initial values correspoinding to each "
                                     "group, i.e. '0.0 1.0 2.0' means 0-4 = 0.0, "
                                     "4-8 = 1.0 and 8-12 = 2.0 ");
  return params;
}

CPPlasticEnergyMaterial::CPPlasticEnergyMaterial(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
  _variable_size(getParam<unsigned int>("variable_size")),
  _mat_prop_state_var(
    getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
  _Q(getParam<Real>("Q")),
  _q(getParam<Real>("q")),
  // _s0(getParam<Real>("s0")),
  _groups(getParam<std::vector<unsigned int>>("groups")),
  _group_values(getParam<std::vector<Real>>("group_values"))
{

}
Real
CPPlasticEnergyMaterial::computeF()
{
  std::vector<Real> s0;
  s0.resize(_variable_size);
  for (unsigned int i = 0; i < _groups.size() - 1; ++i)
  {
    unsigned int is, ie;

    is = _groups[i];
    ie = _groups[i + 1] - 1;

    if (is > ie)
      mooseError("CPPlasticEnergyMaterial: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in state variable read");

    for (unsigned int j = is; j <= ie; ++j)
      s0[j] = _group_values[i];
  }


  Real sum, qab;
  sum = 0;
  // Loop over all slip systems
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < _variable_size; ++j)
    {
      // unsigned int iplane, jplane;
      // iplane = i / 2;
      // jplane = j / 2;
      //
      // // Not sure why not if i==j
      // if (iplane == jplane) // Kalidindi
      if (i == j)
        qab = 1.0;
      // else
      //   qab = _q;
      else if ( i >=20 || j >=20)
        qab = _q;

      else if ((i == 0||i==4)||(i==14||i==15))
        if ((j == 0||j==4)||(j==14||j==15))
          qab = 1.0;
      else if ((i == 1||i==5)||(i==12||i==13))
        if ((j == 1||j==5)||(j==12||j==13))
          qab = 1.0;
      else if ((i==2||i==6)||(i==7||i==10))
        if ((j==2||j==6)||(j==7||j==10))
          qab = 1.0;

      else if ((i==3||i==8)||(i==9||i==11))
        if ((j==3||j==8)||(j==9||j==11))
          qab = 1.0;

      else if (i>=16&&i<=19)
        if (i>=16&&i<=19)
          qab = 1.0;

      else
        qab = _q;

      sum += qab * _mat_prop_state_var[_qp][i] * _mat_prop_state_var[_qp][j];
      // This will not give sum = 0 for purely elastic so must subtract initial value of _mat_prop_state_var
      sum -= qab * s0[i] * s0[j];

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
