#include "KKSMultiACKernelAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"
#include <iostream>
template <>
InputParameters
validParams<KKSMultiACKernelAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription(
      "Sets up TimeDerivative, KKSMultiACBulkF, KKSMultiACBulkC, ACInterface and ACGrGrMulti kernels");
  params.addRequiredParam<unsigned int>(
      "op_num", "specifies the total number of order parameters");
  params.addRequiredParam<std::string>("op_name_base", "specifies the base name of the order parameters");
  params.addRequiredParam<std::string>("ci_name_base", "specifies the base name of the phase concentrations");
  params.addParam<std::string>("f_name_base","fch","specifies the base name of the free energy");
  params.addParam<std::string>("h_name_base","h","specifies the base name of the switching functions");
  params.addParam<std::string>("g_name_base","g","specifies the base name of the barrier functions");
  params.addParam<Real>("wi", 0, "Height of double well");
  params.addParam<MaterialPropertyName>("mob_name","L","Mobility material for AC equations");
  params.addParam<MaterialPropertyName>("kappa","kappa","kappa");
  params.addParam<MaterialPropertyName>("gamma","gamma","gamma (only works for constant gamma)");

  params.addParam<bool>("implicit", true, "Whether kernels are implicit or not");
  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  return params;
}

KKSMultiACKernelAction::KKSMultiACKernelAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _op_name_base(getParam<std::string>("op_name_base")),
    _ci_name_base(getParam<std::string>("ci_name_base")),
    _f_name_base(getParam<std::string>("f_name_base")),
    _h_name_base(getParam<std::string>("h_name_base")),
    _g_name_base(getParam<std::string>("g_name_base")),
    _implicit(getParam<bool>("implicit"))
{
}

void
KKSMultiACKernelAction::act()
{
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    //
    // Create variable names
    //

    std::string var_name = _op_name_base + Moose::stringify(op);
    std::vector<VariableName> eta_i = {var_name}; // "vector" of var_name   {getParam<NonlinearVariableName>("uvel")}
    std::vector<VariableName> v; // vector of the other order parameters
    v.resize(_op_num - 1);
    std::vector<VariableName> args; // vector of the other order parameters and all phase concentrations
    args.resize(2*_op_num - 1);
    std::vector<VariableName> ci; // vector of phase concentrations
    ci.resize(_op_num);
    std::vector<MaterialPropertyName> h; // vector of switching functions
    h.resize(_op_num);
    std::vector<MaterialPropertyName> f; // vector of energies
    f.resize(_op_num);
    std::vector<MaterialPropertyName> gamma; // vector of gamma
    gamma.resize(_op_num-1);

    // std::cout << "Working on variable: "<< var_name <<std::endl;
    // Populate arrays
    unsigned int ind = 0;
    for (unsigned int j = 0; j < _op_num; ++j)
    {
      if (j != op)
      {
        v[ind] = _op_name_base + Moose::stringify(j);
        args[ind] = _op_name_base + Moose::stringify(j);
        args[_op_num-1+ind] = _ci_name_base + Moose::stringify(j);
        gamma[ind] = getParam<MaterialPropertyName>("gamma");
        ind += 1;
      }
      // here we can use j to index because we want to load all
      ci[j] = _ci_name_base + Moose::stringify(j);
      h[j] = _h_name_base + Moose::stringify(j);
      f[j] = _f_name_base + Moose::stringify(j);
    }
    //now we are missing ci_op in ci_j in args so we need to add it
    // args[2*_op_num-2] = _ci_name_base + Moose::stringify(op);
    args[_op_num-1+ind] = _ci_name_base + Moose::stringify(op);


    //
    // Set up TimeDerivative kernels
    //

    {
      InputParameters params = _factory.getValidParams("TimeDerivative");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<bool>("implicit") = true;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");
      std::string kernel_name = "ddt_" + var_name;
      _problem->addKernel("TimeDerivative", kernel_name, params);
    }

    //
    // Set up KKSMultiACBulkF kernel
    //
    {
      InputParameters params = _factory.getValidParams("KKSMultiACBulkF");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("eta_i") = eta_i;
      // params.set<std::vector<VariableName, std::allocator<VariableName>>>("eta_i") = var_name;
      params.set<std::vector<MaterialPropertyName>>("Fj_names") = f;
      params.set<std::vector<MaterialPropertyName>>("hj_names") = h;
      params.set<MaterialPropertyName>("gi_name") = _g_name_base+Moose::stringify(op);
      params.set<Real>("wi") = getParam<Real>("wi");
      params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mob_name");
      params.set<std::vector<VariableName>>("args") = args;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "KKSMultiACBulkF_" + var_name;
      _problem->addKernel("KKSMultiACBulkF", kernel_name, params);
    }
    // //
    // // Set up KKSMultiACBulkC kernel
    // //
    {
      InputParameters params = _factory.getValidParams("KKSMultiACBulkC");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("eta_i") = eta_i;
      params.set<std::vector<MaterialPropertyName>>("Fj_names") = f;
      params.set<std::vector<MaterialPropertyName>>("hj_names") = h;
      params.set<std::vector<VariableName>>("cj_names") = ci;
      params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mob_name");
      params.set<std::vector<VariableName>>("args") = v;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "KKSMultiACBulkC_" + var_name;
      _problem->addKernel("KKSMultiACBulkC", kernel_name, params);
    }

    // Set up ACInterface kernels
    //
    {
      InputParameters params = _factory.getValidParams("ACInterface");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<MaterialPropertyName>("kappa_name") = getParam<MaterialPropertyName>("kappa");
      params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mob_name");
      params.set<std::vector<VariableName>>("args") = v;
      params.set<bool>("implicit") = getParam<bool>("implicit");
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "ACInt_" + var_name;
      _problem->addKernel("ACInterface", kernel_name, params);
    }

    //
    // Set up ACGrGrMulti kernels
    //
    {
      InputParameters params = _factory.getValidParams("ACGrGrMulti");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("v") = v;
      params.set<std::vector<MaterialPropertyName>>("gamma_names") = gamma;
      params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mob_name");
      params.set<std::vector<VariableName>>("args") = v;
      params.set<bool>("implicit") = _implicit;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "ACGrGrMulti_" + var_name;
      _problem->addKernel("ACGrGrMulti", kernel_name, params);
    }
  }
}
