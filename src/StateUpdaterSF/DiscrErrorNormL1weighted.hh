// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium

// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DiscrErrorNormL1weighted_hh
#define ShockFitting_DiscrErrorNormL1weighted_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/StateUpdater.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class define a DiscrErrorNormL1weighted, whose task is to
/// compute the L1 norm of the discretization error weighted on the first 
/// value of the residual over the grid-points of the background mesh

class DiscrErrorNormL1weighted : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  DiscrErrorNormL1weighted(const std::string& objectName);

  /// Destructor
  virtual ~DiscrErrorNormL1weighted();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// compute the norm value
  virtual void update();

private: // functions

  /// assign the variables used in DiscrErrorNormL1weighted to MeshData pattern
  void setMeshData();

  /// assign the variables used in DiscrErrorNormL1weighted to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// value of the residual t=0
  std::vector<double>* firstResidualValue;

  /// Array2D storing the primitive values of the current step
  /// computed in the grid-points of the background mesh
  /// it is assigned to MeshData
  Array2D <double>* primBackgroundMesh;

  /// Array2D storing the primitive values of the previous step
  /// computed in the grid-points of the background mesh
  /// it is assigned to MeshData
  Array2D <double>* primBackgroundMeshOld;   

  /// value of the computed norm
  std::vector<double> normValue;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_DiscrErrorNormL1weighted_hh
