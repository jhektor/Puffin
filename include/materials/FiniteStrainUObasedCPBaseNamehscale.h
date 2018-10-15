/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINUOBASEDCPBASENAMEHSCALE_H
#define FINITESTRAINUOBASEDCPBASENAMEHSCALE_H

// #include "ComputeStressBase.h"
//
// #include "CrystalPlasticitySlipRateBaseName.h"
// #include "CrystalPlasticitySlipResistance.h"
// #include "CrystalPlasticityStateVariable.h"
// #include "CrystalPlasticityStateVarRateComponent.h"
#include "FiniteStrainUObasedCPBaseName.h"

class FiniteStrainUObasedCPBaseNamehscale;

template <>
InputParameters validParams<FiniteStrainUObasedCPBaseNamehscale>();

/**
 * FiniteStrainUObasedCPBaseName uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the
 * material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 *
 * Involves 4 different types of user objects that calculates:
 * State variables - update state variable (derive from CrystalPlasticityStateVariable)
 * State variable evolution compoment - individual component of state variable incremental rate
 * (derive from CrystalPlasticityStateVariableEvolutionRateComponent)
 * Slip resistance - calcuate slip resistances (derive from CrystalPlasticitySlipResistances)
 * Slip rates - calcuate flow direction and slip rates (derive from CrystalPlasticitySlipRates)

 * This is the same code as FiniteStrainUOBasedCP except the base_name option has been added to the following variables:
 * _fp, _fp_old, _lag_e, _pk2, _pl2_old,_lag_e, _update_rot, _update_rot_old, _deformation_gradient, _deformation_gradient_old, _crysrot
 */
class FiniteStrainUObasedCPBaseNamehscale : public FiniteStrainUObasedCPBaseName
{
public:
  FiniteStrainUObasedCPBaseNamehscale(const InputParameters & parameters);

protected:
  /**
   * updates the stress at a quadrature point.
   */
  virtual void computeQpStress();

  // Solves only where this h material is >0
  const std::string _h_scale_name;
  const MaterialProperty<Real> &  _h_scale;
  const Real _h_tol;
};

#endif // FINITESTRAINUOBASEDCPBASENAME_H
