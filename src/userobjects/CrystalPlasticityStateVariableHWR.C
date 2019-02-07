//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityStateVariableHWR.h"

#include <fstream>

registerMooseObject("PuffinApp", CrystalPlasticityStateVariableHWR);

template <>
InputParameters
validParams<CrystalPlasticityStateVariableHWR>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariable>();
  return params;
}

CrystalPlasticityStateVariableHWR::CrystalPlasticityStateVariableHWR(const InputParameters & parameters)
  : CrystalPlasticityStateVariable(parameters)
{
  if (_scale_factor.size() != _num_mat_state_var_evol_rate_comps)
    mooseError("CrystalPlasticityStateVariableHWR: Scale factor should be have the same size of "
               "evolution rate components.");

  _mat_prop_state_var_evol_rate_comps.resize(_num_mat_state_var_evol_rate_comps);

  for (unsigned int i = 0; i < _num_mat_state_var_evol_rate_comps; ++i)
    _mat_prop_state_var_evol_rate_comps[i] = &getMaterialProperty<std::vector<Real>>(
        parameters.get<std::vector<std::string>>("uo_state_var_evol_rate_comp_name")[i]);
}

void
CrystalPlasticityStateVariableHWR::initSlipSysProps(std::vector<Real> & val,
                                                 const Point & q_point) const
{
  switch (_intvar_read_type)
  {
    case 0:
      readInitialValueFromFile(val);
      break;
    case 1:
      readInitialValueFromInline(val);
      break;
    case 2:
      provideInitialValueByUser(val, q_point);
      break;
    default:
      mooseError("CrystalPlasticityStateVariableHWR: Read option for initial value of internal "
                 "variables is not supported.");
  }

  // for (unsigned int i = 0; i < _variable_size; ++i)
  //   if (val[i] <= 0.0)
  //     mooseError("CrystalPlasticityStateVariableHWR: Value of state variables ", i, " non positive");
}
