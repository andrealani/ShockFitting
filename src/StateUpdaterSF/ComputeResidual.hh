// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeResidual_hh
#define ShockFitting_ComputeResidual_hh

//--------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"
#include "Framework/VariableTransformer.hh"
#include "MathTools/Array2D.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeResidual, whose task is to compute
/// the norm L1 of the discretization error using the quantities at the
/// current step and at the previous one

class ComputeResidual : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ComputeResidual(const std::string& objectName);

  /// Destructor
  virtual ~ComputeResidual();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// update phantom nodes values
  virtual void update();

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper functions

  /// assign variables used in ComputeResidual to the MeshData pattern
  void setMeshData();

  /// assign variables used in ComputeResidual to the PhysicsData pattern
  void setPhysicsData();

  /// resize local vector and arrays
  void resizeArray();

  /// assign starting pointer for the array2D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// norm value get from the norm computing object
  double normValue;

  /// number of mesh points
  std::vector <unsigned>* npoin;

  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// working vectors used to exchange data with the object
  /// transforming the variables
  std::vector<double> m_zroe;
  std::vector<double> m_prim;
  std::vector<double> m_XY;

  /// Array2D of the zroe values belonging to the current step
  Array2D <double>* zroe;

  /// Array2D storing the primitive values of the current step
  /// computed in the grid-points of the background mesh
  /// it is assigned to MeshData
  Array2D <double>* primBackgroundMesh;

  /// Array2D storing the primitive values of the previous step
  /// computed in the grid-points of the background mesh
  /// it is assigned to MeshData
  Array2D <double>* primBackgroundMeshOld;

  /// Array2D storing the primitive values of the state in the
  /// grid-points of the background mesh
  /// specifies which kind of norm will be used
  std::string m_whichNorm;

  /// specifies if the used norm is weighted
  bool m_isItWeighted;

  /// specifies the gas model used fot the coconverison in primitive variables
  std::string m_gasModel;

  /// specifies the lowest value of the residual to stop the simulation
  double m_minResidual;

  /// specifies what happens when the lowest residual is reached
  /// by default the "endSimulation" is set
  std::string m_stopAdditionalInfo;

  /// command object computing norm of the discretization error
  PAIR_TYPE(StateUpdater) m_normErr;

  /// command object making the variable conversion. From param to prim dimensional
  PAIR_TYPE(VariableTransformer) m_paramToprimDimensional;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeResidual_hh


