/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSLIPRATEHWR_H
#define CRYSTALPLASTICITYSLIPRATEHWR_H

#include "CrystalPlasticitySlipRateBaseName.h"
#include "RankTwoTensor.h"

class CrystalPlasticitySlipRateHWR;

template <>
InputParameters validParams<CrystalPlasticitySlipRateHWR>();

/**
 * Slip rate userobject class. This class uses the base_name parameter to get _pk2.
 * The model is from P. Håkansson et al. International Journal of Solids and Structures 45 (2008)
 */
class CrystalPlasticitySlipRateHWR : public CrystalPlasticitySlipRateBaseName
{
public:
  CrystalPlasticitySlipRateHWR(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real /*dt*/, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp,
                                 std::vector<RankTwoTensor> & flow_direction) const;

protected:

  const MaterialProperty<std::vector<Real>> & _mat_prop_slip_res;

  const MaterialProperty<RankTwoTensor> & _pk2;

  const Real _gd0; // Initial slip rates
  const Real _m; // rate exponents

  const MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;
};

#endif // CRYSTALPLASTICITYSLIPRATEHWR_H
