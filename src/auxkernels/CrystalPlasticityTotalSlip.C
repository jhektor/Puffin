
#include "CrystalPlasticityTotalSlip.h"
registerMooseObject("PuffinApp", CrystalPlasticityTotalSlip);

template <>
InputParameters
validParams<CrystalPlasticityTotalSlip>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the total slip in the current time step");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("slip_rates", "List of slip rates to sum up. Place in same order as h_names!");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("h_names", "List of switching functions for weighting the slips from different grains. Place in same ofrder as slip_rates");

  return params;
}

CrystalPlasticityTotalSlip::CrystalPlasticityTotalSlip(const InputParameters & parameters)
  : AuxKernel(parameters),
    _h_names(getParam<std::vector<MaterialPropertyName>>("h_names")),
    _slip_rates_names(getParam<std::vector<MaterialPropertyName>>("slip_rates")),
    _num_j(_h_names.size()),
    _prop_h(_num_j),
    _prop_slip_rates(_num_j)
{
  // Check to ensure size of h_names is the same as slip_rates
  if (_slip_rates_names.size() != _num_j)
    mooseError("Size of slip_rates is not equal to the size of h_names in "
               "CrystalPlasticityTotalSlip AuxKernel ",
               name());

  // Get h values
  for (unsigned int i = 0; i < _num_j; ++i)
  {
    _prop_h[i] = &getMaterialPropertyByName<Real>(_h_names[i]);
    _prop_slip_rates[i] = &getMaterialPropertyByName<std::vector<Real>>(_slip_rates_names[i]);
  }
}

Real
CrystalPlasticityTotalSlip::computeValue()
{
  Real total_slip = 0;
  Real slip_nrm;
  Real hi;
  std::vector<Real> slip_rates;
  for (unsigned int i = 0; i < _num_j; ++i)
  {
    slip_rates = (*_prop_slip_rates[i])[_qp];
    slip_nrm = 0;
    for (unsigned int j = 0; j<slip_rates.size(); ++j)
    {
      slip_nrm += slip_rates[j]*slip_rates[j];
    }
    slip_nrm = sqrt(slip_nrm);
    hi = (*_prop_h[i])[_qp];
    total_slip += _dt*hi*slip_nrm;
  }
  return total_slip;
}
