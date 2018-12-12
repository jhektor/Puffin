//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TricrystalTripleJunctionICyscale.h"
// #include "FunctionIC.h"
// #include "TricrystalTripleJunctionIC.h"

registerMooseObject("PuffinApp", TricrystalTripleJunctionICyscale);

template <>
InputParameters
validParams<TricrystalTripleJunctionICyscale>()
{
  InputParameters params = validParams<TricrystalTripleJunctionIC>();
  params.addClassDescription("Tricrystal with a triple junction, allows for triple junction not spanning the full y range");
  // params.addParam<Real>("ylim",-1000,"lowest y coordinate where this IC sets values different from 0");
  params.addRequiredParam<FunctionName>("ylim","function specifying the lowest y coordinate where this IC sets values different from 0");
  return params;
}

TricrystalTripleJunctionICyscale::TricrystalTripleJunctionICyscale(const InputParameters & parameters)
  : TricrystalTripleJunctionIC(parameters),
    // _ylim(getParam<Real>("ylim"))
    _ylim(getFunction("ylim"))
{
}

Real
TricrystalTripleJunctionICyscale::value(const Point & p)
{
  Real icval;
  if (p(1) >= _ylim.value(_t,p))
  {
    icval = TricrystalTripleJunctionIC::value(p);
    return icval; //TricrystalTripleJunctionICyscale::value(p);
  }
  else
    return 0.0;
}
