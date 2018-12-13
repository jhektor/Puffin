#ifndef PUFFINAPP_H
#define PUFFINAPP_H

#include "MooseApp.h"

class PuffinApp;

template <>
InputParameters validParams<PuffinApp>();

class PuffinApp : public MooseApp
{
public:
  PuffinApp(InputParameters parameters);
  virtual ~PuffinApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* PuffinAPP_H  */
