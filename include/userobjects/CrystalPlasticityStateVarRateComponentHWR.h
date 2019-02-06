#ifndef CRYSTALPLASTICITYSTATEVARRATECOMPONENTHWR_H
#define CRYSTALPLASTICITYSTATEVARRATECOMPONENTHWR_H

#include "CrystalPlasticityStateVarRateComponent.h"

class CrystalPlasticityStateVarRateComponentHWR;

template <>
InputParameters validParams<CrystalPlasticityStateVarRateComponentHWR>();

/**
 * Phenomenological constitutive model state variable evolution rate component userobject class.
 */
class CrystalPlasticityStateVarRateComponentHWR : public CrystalPlasticityStateVarRateComponent
{
public:
  CrystalPlasticityStateVarRateComponentHWR(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp,
                                                       std::vector<Real> & val) const;

protected:
  const std::string _base_name;

  const MaterialProperty<std::vector<Real>> & _mat_prop_slip_rate;
  const MaterialProperty<std::vector<Real>> & _mat_prop_state_var;
  const MaterialProperty<std::vector<Real>> & _mat_prop_slip_res;

  const std::vector<Real> _B;

  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;


};

#endif // CRYSTALPLASTICITYSTATEVARRATECOMPONENTHWR_H
