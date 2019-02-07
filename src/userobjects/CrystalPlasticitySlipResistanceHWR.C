//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticitySlipResistanceHWR.h"

registerMooseObject("PuffinApp", CrystalPlasticitySlipResistanceHWR);

template <>
InputParameters
validParams<CrystalPlasticitySlipResistanceHWR>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addParam<std::string>("uo_state_var_name",
                               "Name of state variable property: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<unsigned int>(
                               "crystal_lattice_type",
                               2, // could choose FCC, BCC, or tin BCT crystals
                               "Type of crystal lattice structure: 0: FCC, 1: BCC, 2: BCT_tin");
  params.addParam<std::vector<unsigned int>>("groups",
                                            "To group the initial values on different "
                                            "slip systems 'format: [start end)', i.e.'0 "
                                            "4 8 11' groups 0-3, 4-7 and 8-11 ");
  params.addParam<std::vector<Real>>("G0",
                                    "The initial values of G0 correspoinding to each "
                                    "group, i.e. '0.0 1.0 2.0' means 0-4 = 0.0, "
                                    "4-8 = 1.0 and 8-12 = 2.0 ");
  params.addParam<Real>("Q",1.0,"Material parameter for slip resistance");
  params.addParam<Real>("q",1.0,"ratio between self and latent hardening");
  params.addClassDescription("HWR crystal plasticity model slip resistance class."
                            );
  return params;
}

CrystalPlasticitySlipResistanceHWR::CrystalPlasticitySlipResistanceHWR(
    const InputParameters & parameters)
  : CrystalPlasticitySlipResistance(parameters),
    _mat_prop_state_var(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
    _crystal_lattice_type(getParam<unsigned int>("crystal_lattice_type")),
    _Q(getParam<Real>("Q")),
    _q(getParam<Real>("q")),
    _groups(getParam<std::vector<unsigned int>>("groups")),
    _group_values(getParam<std::vector<Real>>("G0"))
{
  // fill _G0 vector from groups and group_values
  _G0.resize(_groups[_groups.size()-1]);
  if (_groups.size() <= 0)
    mooseError("CrystalPlasticitySlipResistanceHWR: Error in reading initial state variable values: "
               "Specify input in .i file");
  else if (_groups.size() != (_group_values.size() + 1))
    mooseError(
        "CrystalPlasticitySlipResistanceHWR: The size of the groups and group_values does not match.");

  for (unsigned int i = 0; i < _groups.size() - 1; ++i)
  {
    unsigned int is, ie;

    is = _groups[i];
    ie = _groups[i + 1] - 1;

    if (is > ie)
      mooseError("CrystalPlasticitySlipResistanceHWR: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in state variable read");

    for (unsigned int j = is; j <= ie; ++j)
      _G0[j] = _group_values[i];

  }



}

// MooseEnum
// CrystalPlasticitySlipResistanceHWR::crystalLatticeTypeOptions()
// {
//   return MooseEnum("FCC BCC BCT_tin", "BCT_tin");
// }

bool
CrystalPlasticitySlipResistanceHWR::calcSlipResistance(unsigned int qp,
                                                       std::vector<Real> & val) const
{
  Real hab;
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real sum = 0;
    for (unsigned int j = 0; j < _variable_size; ++j)
    {
      switch (_crystal_lattice_type)
      {
        case 0: // FCC
        {
            unsigned int iplane, jplane;
            iplane = i / 3;
            jplane = j / 3;
            if (iplane == jplane) // Kalidindi
              hab = 1.0;
            else
              hab = _q;
            break;
        }
        case 1: // BCC
        {
          mooseError("Slip resistance for BCC is not yet implemented");
        }
        case 2: // BCT_tin
        {
          // This requires a slip system file which is sorted on slip plane normals (5 groups of 4 and 12 singles)
          unsigned int iplane, jplane;
          iplane = i / 4;
          jplane = j / 4;
          if (i==j)
            hab = 1.0;
          else if (iplane == jplane && (i<20 && j<20))
            hab = 1.0;
          else
            hab = _q;
          break;
        }
        default: mooseError("You've requested the wrong lattice type");
      }
      sum += hab*_mat_prop_state_var[qp][j];
    }
    val[i] = _G0[i] + _Q*sum;
  }
  return true;
}
