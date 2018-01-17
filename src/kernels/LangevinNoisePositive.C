#include "LangevinNoisePositive.h"
#include "MooseRandom.h"

template <>
InputParameters
validParams<LangevinNoisePositive>()
{
  InputParameters params = validParams<LangevinNoise>();
  params.addClassDescription("Source term for non-conserved Langevin noise, the noise is always >=0");
  return params;
}
LangevinNoisePositive::LangevinNoisePositive(const InputParameters & parameters)
  : LangevinNoise(parameters)
{

}

Real
LangevinNoisePositive::computeQpResidual()
{
  return -_test[_i][_qp] *  MooseRandom::rand()  * _amplitude * _multiplier_prop[_qp];
}
