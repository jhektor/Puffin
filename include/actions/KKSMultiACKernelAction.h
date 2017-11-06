#ifndef KKSMULTIACKERNELACTION_H
#define KKSMULTIACKERNELACTION_H

#include "Action.h"

/**
 * Action that sets up TimeDerivative, KKSMultiACBulkF, KKSMultiACBulkC, ACInterface and ACGrGrMulti
 * kernels.
 */
class KKSMultiACKernelAction : public Action
{
public:
  KKSMultiACKernelAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _op_name_base;
  const std::string _ci_name_base;
  const std::string _fch_name_base;
  const std::string _h_name_base;
  const std::string _g_name_base;

  /// kernels are implicit?
  const bool _implicit;
};

template <>
InputParameters validParams<KKSMultiACKernelAction>();

#endif // KKSMULTIACKERNELACTION_H
