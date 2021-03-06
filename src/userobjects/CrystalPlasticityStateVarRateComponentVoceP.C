//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityStateVarRateComponentVoceP.h"
#include "MooseError.h"
#include <iostream>  // I/O


registerMooseObject("PuffinApp", CrystalPlasticityStateVarRateComponentVoceP);

template <>
InputParameters
validParams<CrystalPlasticityStateVarRateComponentVoceP>()
{

  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>(
      "uo_slip_rate_name",
      "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addParam<std::string>("uo_state_var_name",
                               "Name of state variable property: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<MooseEnum>(
      "crystal_lattice_type",
      CrystalPlasticityStateVarRateComponentVoceP::crystalLatticeTypeOptions(),
      "Type of crystal lattice structure output");
  params.addParam<std::vector<unsigned int>>("groups",
                                             "To group the initial values on different "
                                             "slip systems 'format: [start end)', i.e.'0 "
                                             "12 24 48' groups 0-11, 12-23 and 24-48 ");
  params.addParam<std::vector<Real>>("h0_group_values",
                                     "h0 hardning constatn for each group "
                                     " i.e. '0.0 1.0 2.0' means 0-11 = 0.0, "
                                     "12-23 = 1.0 and 24-48 = 2.0 ");
  params.addParam<std::vector<Real>>("tau0_group_values",
                                     "The initial critical resolved shear stress"
                                     "correspoinding to each group"
                                     " i.e. '100.0 110.0 120.0' means 0-11 = 100.0, "
                                     "12-23 = 110.0 and 24-48 = 120.0 ");
  params.addParam<std::vector<Real>>("tauSat_group_values",
                                     "The saturation resolved shear stress"
                                     "correspoinding to each group"
                                     " i.e. '150.0 170.0 180.0' means 0-11 = 150.0, "
                                     "12-23 = 170.0 and 24-48 = 180.0 ");
  params.addParam<std::vector<Real>>("hardeningExponent_group_values",
                                     "The hardening exponent m"
                                     "correspoinding to each group"
                                     " i.e. '1.0 2.0 3.0' means 0-11 = 1.0, "
                                     "12-23 = 2.0 and 24-48 = 3.0 ");
  params.addParam<std::vector<Real>>("selfHardening_group_values",
                                     "The self hardening coefficient q_aa"
                                     "correspoinding to each group"
                                     " i.e. '1.0 2.0 3.0' means 0-11 = 1.0, "
                                     "12-23 = 2.0 and 24-48 = 3.0 "
                                     " usually these are all 1.");
  params.addParam<std::vector<Real>>("coplanarHardening_group_values",
                                     "The coplanar laten hardening coefficient q_ab"
                                     "correspoinding to each group"
                                     " i.e. '1.0 2.0 3.0' means 0-11 = 1.0, "
                                     "12-23 = 2.0 and 24-48 = 3.0 ");
  // params.addParam<bool>("_check_for_errors", "Enables error check output in file.");
  // params.addParam<std::vector<Real>>("GroupGroup_Hardening_group_values",
  //                                    "The group-to-group laten hardening coefficient q_ab"
  //                                    "This is a NxN vector"
  //                                    " i.e. '1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0' "
  //                                    "means non-coplanar slip systems in gr_11,22,33= "
  //                                    "1.0, 5.0 and 9.0 respectively."
  //                                    "latent ahrdening between for gr_12,13 = 2.0 3.0"
  //                                    " respectively");
  params.addClassDescription(
      "Phenomenological Voce constitutive model state"
      " variable evolution rate "
      "component base class.  Override this virtual functions in your class");
  return params;
}

CrystalPlasticityStateVarRateComponentVoceP::CrystalPlasticityStateVarRateComponentVoceP(
    const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_slip_rate(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_slip_rate_name"))),
    _mat_prop_state_var(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_state_var_name"))),
    _crystal_lattice_type(getParam<MooseEnum>("crystal_lattice_type")),
    _groups(getParam<std::vector<unsigned int>>("groups")),
    _h0_group_values(getParam<std::vector<Real>>("h0_group_values")),
    _tau0_group_values(getParam<std::vector<Real>>("tau0_group_values")),
    _tauSat_group_values(getParam<std::vector<Real>>("tauSat_group_values")),
    _hardeningExponent_group_values(getParam<std::vector<Real>>("hardeningExponent_group_values")),
    _selfHardening_group_values(getParam<std::vector<Real>>("selfHardening_group_values")),
    _coplanarHardening_group_values(getParam<std::vector<Real>>("coplanarHardening_group_values"))
    // _check_for_errors(false) // change if error check is neccessary
    // _GroupGroup_Hardening_group_values(
    //     getParam<std::vector<Real>>("GroupGroup_Hardening_group_values"))
{
  // perform input checks and initialize usefull variables
  _check_for_errors=false;
  _n_groups = _groups.size() - 1;
  checkHardeningParametersSize();
  initSlipSystem_PlaneID(_slipSystem_PlaneID);
  initSlipSystem_GroupID(_slipSystem_GroupID);
}

void
CrystalPlasticityStateVarRateComponentVoceP::checkHardeningParametersSize() const
{
  // check that at least one group exists
  if (_n_groups <= 0)
    mooseError("CrystalPlasticityStateVarRateComponentVoce: Error in reading hardening     "
               "parameters values "
               "Specify input in .i file ");

  // check the size of all the variables
  bool check_var_size = true;
  check_var_size &= _h0_group_values.size() == _n_groups;
  check_var_size &= _tau0_group_values.size() == _n_groups;
  check_var_size &= _tauSat_group_values.size() == _n_groups;
  check_var_size &= _hardeningExponent_group_values.size() == _n_groups;
  check_var_size &= _selfHardening_group_values.size() == _n_groups;
  check_var_size &= _coplanarHardening_group_values.size() == _n_groups;
  check_var_size &= _h0_group_values.size() == _n_groups;
  // check_var_size &= _GroupGroup_Hardening_group_values.size() == _n_groups * _n_groups;

  if (!check_var_size)
    mooseError("CrystalPlasticityStateVarRateComponentVoce: "
               "The size of one or more input parameters does not match the group size");
}

bool
CrystalPlasticityStateVarRateComponentVoceP::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  unsigned int group_i;
  Real h0;
  Real tau_0;
  Real tau_sat;
  Real hardening_exponenet;
  Real delta_tau;

  DenseVector<Real> hb(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    group_i = _slipSystem_GroupID[i];
    h0 = _h0_group_values[group_i];
    tau_0 = _tau0_group_values[group_i];
    tau_sat = _tauSat_group_values[group_i];
    hardening_exponenet = _hardeningExponent_group_values[group_i];

    delta_tau = tau_sat - tau_0;

    hb(i) = h0 *
            std::pow(std::abs(1.0 - (_mat_prop_state_var[qp][i] - tau_0) / delta_tau),
                     hardening_exponenet) *
            std::copysign(1.0, 1.0 - (_mat_prop_state_var[qp][i] - tau_0) / delta_tau);
  }

  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < _variable_size; ++j)
    {

      Real q_ab = getHardeningCoefficient(i, j);

      val[i] += std::abs(_mat_prop_slip_rate[qp][j]) * q_ab * hb(j);

    }

  return true;
}

MooseEnum
CrystalPlasticityStateVarRateComponentVoceP::crystalLatticeTypeOptions()
{
  return MooseEnum("FCC BCC BCT", "BCT"); // Inserted BCT into list and as default 8/5-18
}

void
CrystalPlasticityStateVarRateComponentVoceP::initSlipSystem_PlaneID(
    std::vector<unsigned int> & _slipSystem_PlaneID) const
{
  // this routine is generate a vector containing the association between
  // slip system number and slip plane
  _slipSystem_PlaneID.assign(_variable_size, 0);
  if (_check_for_errors)
    std::cerr << " # Running: 'initSlipSystem_PlaneID' # _variable_size = " << _variable_size << std::endl;

  for (unsigned int slipSystemIndex = 0; slipSystemIndex < _variable_size; ++slipSystemIndex)
  {

    switch (_crystal_lattice_type)
    {
      case 0: // FCC
        if (slipSystemIndex < 12)
          _slipSystem_PlaneID[slipSystemIndex] = slipSystemIndex / 3;
        else
          mooseError("FCC with more than 12 slip planes not implemented ");

        break;

      case 1: // BCC
        if (slipSystemIndex < 12)
          _slipSystem_PlaneID[slipSystemIndex] = slipSystemIndex / 2;

        else if (slipSystemIndex >= 12 && slipSystemIndex < 48)
          _slipSystem_PlaneID[slipSystemIndex] = (slipSystemIndex - 6);

        else
          mooseError("BCC with more than 48 slip systems not implemented ");

        break;
      case 2: // BCT, see 'slip_systems_bct_v2_sorted2.txt' for details
        // if (slipSystemIndex < 2)
        //   _slipSystem_PlaneID[slipSystemIndex] = 0;
        //
        // else if (slipSystemIndex < 4)
        //   _slipSystem_PlaneID[slipSystemIndex] = 1;
        //
        // else if (slipSystemIndex < 6)
        //   _slipSystem_PlaneID[slipSystemIndex] = 2;
        //
        // else if (slipSystemIndex < 10)
        //   _slipSystem_PlaneID[slipSystemIndex] = 3;
        //
        // else if (slipSystemIndex < 20)
        //   _slipSystem_PlaneID[slipSystemIndex] = 4;
        //
        // else if (slipSystemIndex < 24)
        //   _slipSystem_PlaneID[slipSystemIndex] = 5;
        //
        // else if (slipSystemIndex >= 24 && slipSystemIndex < 32)
        //   _slipSystem_PlaneID[slipSystemIndex] = 6;
        if (slipSystemIndex >=20 )
          _slipSystem_PlaneID[slipSystemIndex] = slipSystemIndex;

        else if ((slipSystemIndex == 0||slipSystemIndex==4)||(slipSystemIndex==14||slipSystemIndex==15))
          _slipSystem_PlaneID[slipSystemIndex] = 0;

        else if ((slipSystemIndex == 1||slipSystemIndex==5)||(slipSystemIndex==12||slipSystemIndex==13))
          _slipSystem_PlaneID[slipSystemIndex] = 1;

        else if ((slipSystemIndex==2||slipSystemIndex==6)||(slipSystemIndex==7||slipSystemIndex==10))
          _slipSystem_PlaneID[slipSystemIndex] = 2;

        else if ((slipSystemIndex==3||slipSystemIndex==8)||(slipSystemIndex==9||slipSystemIndex==11))
          _slipSystem_PlaneID[slipSystemIndex] = 3;

        else if (slipSystemIndex>=16&&slipSystemIndex<=19)
          _slipSystem_PlaneID[slipSystemIndex] = 4;

        else
          mooseError("There is something shady in the system Sire! BCT slip systems are compromised");

        if (_check_for_errors) //only for debugging
          std::cerr << "# value of slipSystemIndex " << slipSystemIndex << " # planeID of ss " << _slipSystem_PlaneID[slipSystemIndex]<< std::endl;

        break;
      default:
        mooseError("VoceHardeningError: Pass valid crustal_structure_type ");
    }
  }
}

void
CrystalPlasticityStateVarRateComponentVoceP::initSlipSystem_GroupID(
    std::vector<unsigned int> & _slipSystem_GroupID) const
// this routine is generate a vector containing the association between
// slip system number and provided group edges
{
  _slipSystem_GroupID.assign(_variable_size, 0);

  for (unsigned int slipSystemIndex = 0; slipSystemIndex < _variable_size; ++slipSystemIndex)
  {
    for (unsigned int i = 0; i < _groups.size() - 1; i++)
    {
      if (slipSystemIndex >= _groups[i] && slipSystemIndex < _groups[i + 1])
      {
        _slipSystem_GroupID[slipSystemIndex] = i;
        break;
      }
    }
  }
}

Real
CrystalPlasticityStateVarRateComponentVoceP::getHardeningCoefficient(
    unsigned int slipSystemIndex_i, unsigned int slipSystemIndex_j) const
{
  // select the appropriate latent hardening coefficient based on the slip systems indices

  Real q_ab;
  // collect slip system plane and group
  unsigned int group_i = _slipSystem_GroupID[slipSystemIndex_i];
  unsigned int group_j = _slipSystem_GroupID[slipSystemIndex_j];
  unsigned int plane_i = _slipSystem_PlaneID[slipSystemIndex_i];
  unsigned int plane_j = _slipSystem_PlaneID[slipSystemIndex_j];

  // create check for clarity
  bool same_slipSystem = slipSystemIndex_i == slipSystemIndex_j;
  bool same_group = group_i == group_j;
  bool same_plane = plane_i == plane_j;

  // retrieve appropriate coefficient
  if (same_slipSystem)
    q_ab = _coplanarHardening_group_values[group_i];
  else if (same_plane)
    q_ab = _coplanarHardening_group_values[group_i];
  // else // here for debugging purposes
  //   mooseError("VoceHardeningError:getHardeningCoefficient: case not listed, abort ");
  else
    q_ab = _selfHardening_group_values[group_i]; //  make this less hacky

  return q_ab;
}
