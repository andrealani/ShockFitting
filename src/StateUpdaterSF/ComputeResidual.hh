// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeResidual_hh
#define ShockFitting_ComputeResidual_hh

//--------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"
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

  /// assign starting pointer for the array2D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// backup of the old number of shocked mesh points
  unsigned* npoinShockedMeshBkp;

  /// norm value get from the norm computing object
  double normValue;

  /// number of mesh points
  std::vector <unsigned>* npoin;

  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points state belonging to the previous time step (assignable to MeshData)
  std::vector<double>* zroeOldVect;

  /// Array2D of the zroe values belonging to thge current step
  Array2D <double>* zroe;

  /// Array2D of the zroe values belonging to the previous step
  Array2D <double>* zroeOld;

  /// specifies which kind of norm will be used
  std::string m_whichNorm;

  /// specifies if the used norm is weighted
  bool m_isItWeighted;

  /// command object computing norm of the discretization error
  PAIR_TYPE(StateUpdater) m_normErr;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeResidual_hh


