#ifndef LANGEVINNOISEPOSITIVE_H
#define LANGEVINNOISEPOSITIVE_H

// #include "Kernel.h"
#include "LangevinNoise.h"

// Forward Declarations
class LangevinNoisePositive;

template <>
InputParameters validParams<LangevinNoisePositive>();

class LangevinNoisePositive : public LangevinNoise
{
public:
  LangevinNoisePositive(const InputParameters & parameters);

protected:
//  virtual void residualSetup();
  virtual Real computeQpResidual() override;

  // const Real _amplitude;
  // const MaterialProperty<Real> & _multiplier_prop;
};

#endif // LANGEVINNOISEPOSITIVE_H
