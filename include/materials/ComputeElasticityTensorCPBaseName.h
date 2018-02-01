/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEELASTICITYTENSORCPBASENAME_H
#define COMPUTEELASTICITYTENSORCPBASENAME_H

#include "ComputeElasticityTensor.h"
#include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

/**
 * ComputeElasticityTensorCPBaseName defines an elasticity tensor material object for crystal plasticity. This version adds the base_name option to _crysrot
 */
class ComputeElasticityTensorCPBaseName : public ComputeElasticityTensor
{
public:
  ComputeElasticityTensorCPBaseName(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  virtual void assignEulerAngles();

  /**
   * Element property read user object
   * Presently used to read Euler angles -  see test
   */
  const ElementPropertyReadFile * _read_prop_user_object;

  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;

  /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot;

  /// Rotation matrix
  RotationTensor _R;

  // const std::string _base_name;

};

#endif // COMPUTEELASTICITYTENSORCPBASENAME_H
