#ifndef CRYSTALPLASTICITYTOTALSLIP_H
#define CRYSTALPLASTICITYTOTALSLIP_H

#include "AuxKernel.h"
#include "Material.h"

// Forward Declarations
class CrystalPlasticityTotalSlip;

template <>
InputParameters validParams<CrystalPlasticityTotalSlip>();

/**
 * Compute the total slip in the current timestep
 * $\gamma^{tot}=\sum_i h_i\Delta t\sqrt{\sum_k \dot{\gamma}^k_i}\dot{\gamma}^k_i}
 */
class CrystalPlasticityTotalSlip : public AuxKernel
{
public:
  CrystalPlasticityTotalSlip(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// Switching function names
  std::vector<MaterialPropertyName> _h_names;

  /// Slip rates  names
  std::vector<MaterialPropertyName> _slip_rates_names;

  /// Number of grains to interpolate between
  const unsigned int _num_j;

  /// Values of the switching functions for each phase \f$ h_j \f$
  std::vector<const MaterialProperty<Real> *> _prop_h;

  /// Values of the slip_rates
  std::vector<const MaterialProperty<std::vector<Real>> *> _prop_slip_rates;



};

#endif // CRYSTALPLASTICITYTOTALSLIP_H
