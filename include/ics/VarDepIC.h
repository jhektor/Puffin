/**
This initial condition set the varlue of the variable c as c=\sum c_i*eta_i where c_i and eta_i are inputs. The header is modified from MTICSum.h
**/
#ifndef VARDEPIC_H
#define VARDEPIC_H

#include "InitialCondition.h"

class VarDepIC;

template <>
InputParameters validParams<VarDepIC>();

/**
 *
 */
class VarDepIC : public InitialCondition
{
public:
  VarDepIC(const InputParameters & parameters);
  virtual ~VarDepIC();

  virtual Real value(const Point & /*p*/);

protected:
  unsigned int _num_eta; // number of order parameters
  std::vector<const VariableValue *> _etas; // order parameter values
  std::vector<const VariableValue *> _cis; // phase concentration values
};

#endif /* VarDepIC_H */
