#include "UnitySubVarIC.h"

template <>
InputParameters
validParams<UnitySubVarIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredCoupledVar("etas", "Vector of order parameters");
  params.addParam<Real>("y_threshold", 0.0, "Sets variable to 0 at y coordinates below this value");
  params.addParam<bool>("use_threshold", false, "Turns on/off y_threshold");

  return params;
}

UnitySubVarIC::UnitySubVarIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _num_eta(coupledComponents("etas")),
    _etas(_num_eta),
    _y_threshold(parameters.get<Real>("y_threshold")),
    _threshold(parameters.get<bool>("use_threshold"))
{
  // Fetch eta values
  for (unsigned int i = 0; i < _num_eta; ++i)
    _etas[i] = &coupledValue("etas", i);

}

UnitySubVarIC::~UnitySubVarIC() {}

Real
UnitySubVarIC::value(const Point & p)
{
  // If quadrature point is below y_threshold
  if (_threshold && p(1)<_y_threshold)
    return 0.0;

  Real sum_ec = 0.0;
  for (unsigned int i = 0; i < _num_eta; ++i)
    sum_ec += (*_etas[i])[_qp];

  return 1.0-sum_ec;
}
