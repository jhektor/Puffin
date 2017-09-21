/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef IMCFRACTION_H
#define IMCFRACTION_H

#include "NodalVariablePostprocessor.h"

// Forward Declarations
class IMCFraction;

template <>
InputParameters validParams<IMCFraction>();

/**
 * Computes the
 */
class IMCFraction : public NodalVariablePostprocessor
{
public:
  IMCFraction(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

  void threadJoin(const UserObject & y) override;

protected:
  Real _fraction;
  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  std::vector<unsigned int> _vals_var;
};

#endif // IMCFRACTION_H
