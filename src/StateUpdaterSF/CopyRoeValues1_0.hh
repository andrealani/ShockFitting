// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CopyRoeValues1_0_hh
#define ShockFitting_CopyRoeValues1_0_hh

//--------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CopyRoeValues1_0, whose task is to overwrite the 
/// new grid values (1) with the old ones (0).
/// This copy is necessary here since FixStateSps and ComputeStateDps work
/// on nodal values of grid(0)

class CopyRoeValues1_0 : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CopyRoeValues1_0(const std::string& objectName);

  /// Destructor
  virtual ~CopyRoeValues1_0();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// copy old zroe values in the new ones
  virtual void update();

private: // helper functions

  /// assign variables used in CopyRoevalues1_0 to MeshData pattern
  void setMeshData();

  /// assign variables used in CopyRoevalues1_0 to PhysicsData pattern 
  void setPhysicsData();

  /// assign start pointers for arrays 2D
  void setAddress();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// vector of mesh points status (assignable to MeshData pattern)
  std::vector<double>* zroeVect;

  /// array storing new zroe values
  Array2D <double>* zroe1;

  /// array storing old zroe values
  Array2D <double>* zroe0;

  /// map vector
  std::vector <unsigned>* M02M1;

  /// dummy variable to store starting pointer of Array2D
  unsigned start;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_CopyRoeValues1_0_hh

