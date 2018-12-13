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
  params.addParam<MooseEnum>(
                               "crystal_lattice_type",
                               CrystalPlasticitySlipResistanceHWR::crystalLatticeTypeOptions(), // could choose FCC, BCC, or tin BCT crystals
                               "Type of crystal lattice structure output");
  params.addParam<Real>("Q",1.0,"Material parameter for slip resistance");
  params.addParam<Real>("G0",1.0,"Lattice friction part of slip resistance");
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
    _crystal_lattice_type(getParam<MooseEnum>("crystal_lattice_type")),
    _Q(getParam<Real>("Q")),
    _G0(getParam<Real>("G0")),
    _q(getParam<Real>("q"))
{
}

MooseEnum
CrystalPlasticitySlipResistanceHWR::crystalLatticeTypeOptions()
{
  return MooseEnum("FCC BCC BCT_tin", "BCT_tin");
}

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
        }
        case 1: // BCC
          mooseError("Slip resistance for BCC is not yet implemented");
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
        }
      }
      sum += hab*_mat_prop_state_var[qp][j];
    }
    val[i] = _G0 + _Q*sum;
  }
  return true;
}
