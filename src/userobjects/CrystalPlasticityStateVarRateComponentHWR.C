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
     params.addParam<std::vector<unsigned int>>("groups",
                                               "To group the initial values on different "
                                               "slip systems 'format: [start end)', i.e.'0 "
                                               "4 8 11' groups 0-3, 4-7 and 8-11 ");
     params.addParam<std::vector<Real>>("B",
                                       "The initial values of B correspoinding to each "
                                       "group, i.e. '0.0 1.0 2.0' means 0-4 = 0.0, "
                                       "4-8 = 1.0 and 8-12 = 2.0 ");
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
    _groups(getParam<std::vector<unsigned int>>("groups")),
    _group_values(getParam<std::vector<Real>>("B")),
    _pk2(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(parameters.get<std::string>("uo_slip_rate_name")+"_flow_direction"))
{
  // fill _B vector from groups and group_values
  _B.resize(_groups[_groups.size()-1]);
  if (_groups.size() <= 0)
    mooseError("CrystalPlasticityStateVarRateComponentHWR: Error in reading initial state variable values: "
               "Specify input in .i file");
  else if (_groups.size() != (_group_values.size() + 1))
    mooseError(
        "CrystalPlasticityStateVarRateComponentHWR: The size of the groups and group_values does not match.");

  for (unsigned int i = 0; i < _groups.size() - 1; ++i)
  {
    unsigned int is, ie;

    is = _groups[i];
    ie = _groups[i + 1] - 1;

    if (is > ie)
      mooseError("CrystalPlasticityStateVarRateComponentHWR: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in state variable read");

    for (unsigned int j = is; j <= ie; ++j)
      {
        _B[j] = _group_values[i];
      }
  }
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
