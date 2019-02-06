/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CrystalPlasticitySlipRateHWR.h"

#include <fstream>

registerMooseObject("PuffinApp", CrystalPlasticitySlipRateHWR);

template <>
InputParameters
validParams<CrystalPlasticitySlipRateHWR>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRateBaseName>(); // all parameters from here as well.
  params.addParam<std::string>("uo_slip_res_name",
                               "Name of slip resistance property: Same as "
                               "slip resistance user object specified in input "
                               "file."); // plus this one
  params.addParam<Real>("gamma0",1.0, "Initial slip rate");
  params.addParam<Real>("m",1.0, "Rate exponent");
  params.addClassDescription("HWR crystal plastcity model slip rate class.");
  return params;
}

CrystalPlasticitySlipRateHWR::CrystalPlasticitySlipRateHWR(const InputParameters & parameters)
  : CrystalPlasticitySlipRateBaseName(parameters),
    _mat_prop_slip_res(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_slip_res_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk2")),
    _gd0(getParam<Real>("gamma0")),
    _m(getParam<Real>("m")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(_name + "_flow_direction"))
{
}

void
CrystalPlasticitySlipRateHWR::calcFlowDirection(unsigned int qp,
                                                std::vector<RankTwoTensor> & flow_direction) const
{
  DenseVector<Real> mo(LIBMESH_DIM * _variable_size), no(LIBMESH_DIM * _variable_size);

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i * LIBMESH_DIM + j) =
            mo(i * LIBMESH_DIM + j) + _crysrot[qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i * LIBMESH_DIM + j) =
            no(i * LIBMESH_DIM + j) + _crysrot[qp](j, k) * _no(i * LIBMESH_DIM + k);
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        flow_direction[i](j, k) = mo(i * LIBMESH_DIM + j) * no(i * LIBMESH_DIM + k);
}

bool
CrystalPlasticitySlipRateHWR::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  Real tau;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    tau = _pk2[qp].doubleContraction(_flow_direction[qp][i]); // Here it is assumed that the elastic deformation is small

    val[i] = _gd0 * std::pow(std::abs(tau / _mat_prop_slip_res[qp][i]), _m) *
             copysign(1.0, tau);
    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(val[i]) * dt);
#endif
      return false;
    }
  }

  return true;
}

bool
CrystalPlasticitySlipRateHWR::calcSlipRateDerivative(unsigned int qp,
                                                     Real /*dt*/,
                                                     std::vector<Real> & val) const
{
  Real tau;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    tau = _pk2[qp].doubleContraction(_flow_direction[qp][i]);
    val[i] = _gd0 * _m *
             std::pow(std::abs(tau / _mat_prop_slip_res[qp][i]), _m - 1.0) /
             _mat_prop_slip_res[qp][i];
  }
  return true;
}
