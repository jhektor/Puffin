#include "PuffinApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template<>
InputParameters validParams<PuffinApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

PuffinApp::PuffinApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  PuffinApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  PuffinApp::associateSyntax(_syntax, _action_factory);
}

PuffinApp::~PuffinApp()
{
}

// External entry point for dynamic application loading
extern "C" void PuffinApp__registerApps() { PuffinApp::registerApps(); }
void
PuffinApp::registerApps()
{
  registerApp(PuffinApp);
}

// External entry point for dynamic object registration
extern "C" void PuffinApp__registerObjects(Factory & factory) { PuffinApp::registerObjects(factory); }
void
PuffinApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void PuffinApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { PuffinApp::associateSyntax(syntax, action_factory); }
void
PuffinApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
