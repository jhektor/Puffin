//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CrystalPlasticityStateVariableHWR_H
#define CrystalPlasticityStateVariableHWR_H

#include "CrystalPlasticityStateVariable.h"

class CrystalPlasticityStateVariableHWR;

template <>
InputParameters validParams<CrystalPlasticityStateVariableHWR>();

/**
 * Crystal plasticity state variable userobject class.
 */
class CrystalPlasticityStateVariableHWR : public CrystalPlasticityStateVariable
{
public:
  CrystalPlasticityStateVariableHWR(const InputParameters & parameters);
//
//   virtual bool updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual void initSlipSysProps(std::vector<Real> & val, const Point & q_point) const;

};

#endif // CrystalPlasticityStateVariableHWR_H
