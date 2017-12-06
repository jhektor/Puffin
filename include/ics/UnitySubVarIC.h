/**
This initial condition set the value of the variable c as c=1-\sum eta_i where eta_i are inputs. The header is modified from MTICSum.h
**/
#ifndef UNITYSUBVARIC_H
#define UNITYSUBVARIC_H

#include "InitialCondition.h"

class UnitySubVarIC;

template <>
InputParameters validParams<UnitySubVarIC>();

/**
 *
 */
class UnitySubVarIC : public InitialCondition
{
public:
  UnitySubVarIC(const InputParameters & parameters);
  virtual ~UnitySubVarIC();

  virtual Real value(const Point & p);

protected:
  unsigned int _num_eta; // number of order parameters
  std::vector<const VariableValue *> _etas; // order parameter values
  Real _y_threshold;
  bool _threshold;
};

#endif /* UnitySubVarIC_H */
