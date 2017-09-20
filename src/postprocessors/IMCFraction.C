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

#include "IMCFraction.h"
#include "MooseMesh.h"
#include "SubProblem.h"

template <>
InputParameters
validParams<IMCFraction>()
{
  InputParameters params = validParams<NodalVariablePostprocessor>();
  params.set<bool>("unique_node_execute") = true;
  params.addRequiredCoupledVar("eta","Array containing the rest of the order parameters");
  return params;
}

IMCFraction::IMCFraction(const InputParameters & parameters)
  : NodalVariablePostprocessor(parameters),
   _fraction(0),
   _op_num(coupledComponents("eta")), // number of other order parameters
   _vals(_op_num), // array storing the other order parameters
   _vals_var(_op_num) // array with values of the other order parmeters
{
  // Loop through grains and load coupled variables into the arrays
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("eta", i);
    _vals_var[i] = coupled("eta", i);
  }
}

void
IMCFraction::initialize()
{
  _fraction = 0;
}

void
IMCFraction::execute()
{
  // Sum all order parameters
  Real SumEtaj = _u[_qp];
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtaj += (*_vals[i])[_qp];

  _fraction += _u[_qp]/SumEtaj;
}

Real
IMCFraction::getValue()
{
  gatherSum(_fraction);

  return _fraction/_subproblem.mesh().nNodes();
}

void
IMCFraction::threadJoin(const UserObject & y)
{
  const IMCFraction & pps = static_cast<const IMCFraction &>(y);
  _fraction += pps._fraction;
}
