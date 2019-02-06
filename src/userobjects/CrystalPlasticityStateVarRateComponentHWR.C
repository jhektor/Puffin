//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityStateVarRateComponentHWR.h"
#include <cmath>

registerMooseObject("PuffinApp", CrystalPlasticityStateVarRateComponentHWR);

template <>
InputParameters
validParams<CrystalPlasticityStateVarRateComponentHWR>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>(
      "uo_slip_rate_name",
      "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addParam<std::string>("uo_state_var_name",
                               "Name of state variable property: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<std::string>("uo_slip_resistance_name",
                               "Name of slip resistance property: Same as "
                               "slip resistance user object specified in input "
                               "file.");
  params.addRequiredParam<std::vector<Real>>("B", "Material property for saturation of g");
  params.addClassDescription("HWR crystal CrystalPlasticityUOBase model state variable evolution rate ");
  return params;
}

CrystalPlasticityStateVarRateComponentHWR::CrystalPlasticityStateVarRateComponentHWR(
    const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mat_prop_slip_rate(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_slip_rate_name"))),
    _mat_prop_state_var(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
    _mat_prop_slip_res(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_slip_resistance_name"))),
    _B(getParam<std::vector<Real>>("B")),
    _pk2(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(parameters.get<std::string>("uo_slip_rate_name")+"_flow_direction"))
{
  if (_B.size() != _variable_size)
    mooseError("CrystalPlasticityStateVarRateComponentHWR: Size of B does not match the number of slip systems.");
}

bool
CrystalPlasticityStateVarRateComponentHWR::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  Real tau;
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    tau = _pk2[qp].doubleContraction(_flow_direction[qp][i]); // not sure this is correct
    val[i] = (1.0 - _B[i]*_mat_prop_state_var[qp][i])*(std::abs(tau)/_mat_prop_slip_res[qp][i])*std::abs(_mat_prop_slip_rate[qp][i]);
  }
  return true;
}
