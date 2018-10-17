#include "PuffinApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

//PUFFIN
//POSTPROCESSORS
#include "IMCFraction.h"

//InitialConditions
#include "VarDepIC.h"
#include "UnitySubVarIC.h"
#include "TricrystalTripleJunctionICyscale.h"

//Actions
#include "KKSMultiACKernelAction.h"

//Kernels
#include "LangevinNoisePositive.h"

//Materials
#include "ComputeElasticityTensorCPBaseName.h"
#include "FiniteStrainUObasedCPBaseName.h"
#include "FiniteStrainUObasedCPBaseNamehscale.h"
#include "ElasticEnergyMaterialGreenPK2.h"
#include "CPPlasticEnergyMaterial.h"

//UserObjects
#include "CrystalPlasticitySlipRateGSSBaseName.h"
#include "CrystalPlasticitySlipRateBaseName.h"
#include "CrystalPlasticityStateVarRateComponentVoce.h"
// #include "CrystalPlasticitySlipResistanceGSSBaseName.h"
// #include "CrystalPlasticitySlipResistanceBaseName.h"

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
  // Register new stuff here
  registerPostprocessor(IMCFraction);
  registerInitialCondition(VarDepIC);
  registerInitialCondition(UnitySubVarIC);
  registerInitialCondition(TricrystalTripleJunctionICyscale);
  registerKernel(LangevinNoisePositive);
  registerMaterial(ComputeElasticityTensorCPBaseName);
  registerMaterial(FiniteStrainUObasedCPBaseName);
  registerMaterial(FiniteStrainUObasedCPBaseNamehscale);
  registerMaterial(ElasticEnergyMaterialGreenPK2);
  registerMaterial(CPPlasticEnergyMaterial);
  registerUserObject(CrystalPlasticitySlipRateGSSBaseName);
  // registerUserObject(CrystalPlasticityStateVarRateComponentVoce);
  // registerUserObject(CrystalPlasticitySlipResistanceGSSBaseName);

}

// External entry point for dynamic syntax association
extern "C" void PuffinApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { PuffinApp::associateSyntax(syntax, action_factory); }
void
PuffinApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  // Register Action stuff here
  registerAction(KKSMultiACKernelAction, "add_kernel");
  registerSyntax("KKSMultiACKernelAction","Kernels/KKSMultiACKernel");
}
