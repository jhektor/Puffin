#include "PuffinApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<PuffinApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

PuffinApp::PuffinApp(InputParameters parameters) : MooseApp(parameters)
{
  PuffinApp::registerAll(_factory, _action_factory, _syntax);
}

PuffinApp::~PuffinApp() {}

void
PuffinApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"PuffinApp"});
  Registry::registerActionsTo(af, {"PuffinApp"});

  /* register custom execute flags, action syntax, etc. here */
  registerSyntax("KKSMultiACKernelAction","Kernels/KKSMultiACKernel");
}

void
PuffinApp::registerApps()
{
  registerApp(PuffinApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
PuffinApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PuffinApp::registerAll(f, af, s);
}
extern "C" void
PuffinApp__registerApps()
{
  PuffinApp::registerApps();
}
