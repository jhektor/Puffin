//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef TRICRYSTALTRIPLEJUNCTIONICYSCALE_H
#define TRICRYSTALTRIPLEJUNCTIONICYSCALE_H

#include "TricrystalTripleJunctionIC.h"
#include "Function.h"

// Forward Declarations
class TricrystalTripleJunctionICyscale;

template <>
InputParameters validParams<TricrystalTripleJunctionICyscale>();

/**
 * TricrystalTripleJunctionIC creates a 3-grain structure with a triple junction
 * centered at _junction as specified by the user.
 * The initial condition only modifies the values at y > ylim
 */
class TricrystalTripleJunctionICyscale : public TricrystalTripleJunctionIC
{
public:
  TricrystalTripleJunctionICyscale(const InputParameters & parameters);

  virtual Real value(const Point & p);
protected:
  // const Real _ylim;
  Function & _ylim;

};

#endif // TRICRYSTALTRIPLEJUNCTIONICYSCALE_H
