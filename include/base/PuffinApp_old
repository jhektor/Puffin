#ifndef PUFFINAPP_H
#define PUFFINAPP_H

#include "MooseApp.h"

class PuffinApp;

template<>
InputParameters validParams<PuffinApp>();

class PuffinApp : public MooseApp
{
public:
  PuffinApp(InputParameters parameters);
  virtual ~PuffinApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* PUFFINAPP_H */
