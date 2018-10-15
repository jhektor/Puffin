/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainUObasedCPBaseNamehscale.h"
#include "FiniteStrainUObasedCPBaseName.h"
// #include "petscblaslapack.h"
// #include "MooseException.h"
// #include "CrystalPlasticitySlipRateBaseName.h"
// #include "CrystalPlasticitySlipResistance.h"
// #include "CrystalPlasticityStateVariable.h"
// #include "CrystalPlasticityStateVarRateComponent.h"

template <>
InputParameters
validParams<FiniteStrainUObasedCPBaseNamehscale>()
{
  InputParameters params = validParams<FiniteStrainUObasedCPBaseName>();
  params.addParam<std::string>("h_scale","Only solve where h_scale>h_tol");
  params.addParam<Real>("h_tol",0.01,"Only solve where h_scale>h_tol");
  // params
  return params;
}

FiniteStrainUObasedCPBaseNamehscale::FiniteStrainUObasedCPBaseNamehscale(const InputParameters & parameters)
  : FiniteStrainUObasedCPBaseName(parameters),
    _h_scale_name(isParamValid("h_scale") ? getParam<std::string>("h_scale") : ""),
    _h_scale(getMaterialPropertyByName<Real>(_h_scale_name)),
    _h_tol(getParam<Real>("h_tol"))
  {
  }

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
void
FiniteStrainUObasedCPBaseNamehscale::computeQpStress()
{
  // Only update stress if h_scale[_qp]>=h_tol
  if (_h_scale[_qp] <= _h_tol )
  {
    // reset everything
    FiniteStrainUObasedCPBaseName::initQpStatefulProperties();

    // FiniteStrainUObasedCPBaseName::preSolveQp();
    // _dfgrd_tmp = _deformation_gradient[_qp];
    _pk2[_qp].zero();
    // _fe.zero();
    // _fe.addIa(1.0);
    _fp[_qp].zero();
    _fp[_qp].addIa(1.0);
    _dfgrd_tmp.zero();
    _dfgrd_tmp.addIa(1.0);



  }
  else
  {
    FiniteStrainUObasedCPBaseName::computeQpStress();
  }

  FiniteStrainUObasedCPBaseName::postSolveQp();
}
