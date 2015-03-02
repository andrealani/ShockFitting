// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CopyRoeValues2_hh
#define ShockFitting_CopyRoeValues2_hh

//--------------------------------------------------------------------------//

#include "Framework/CopyMaker.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CopyRoeValues2, whose task is to update the 
/// nodal values of the shocks on the new grid 
/// In the fortran version the new grid is referred to index "1"
/// here it is pushed back at the end of the background grid.
/// This copy is necessary here since FixStateSps and ComputeStateDps work
/// on nodal values of grid(0).

class CopyRoeValues2 : public CopyMaker {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CopyRoeValues2(const std::string& objectName);

  /// Destructor
  virtual ~CopyRoeValues2();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// copy old zroe values in the new ones
  virtual void copy();

private: // helper functions

  /// assign variables used in CopyRoevalues1_0 to MeshData pattern
  void setMeshData();

  /// assign variables used in CopyRoevalues1_0 to PhysicsData pattern 
  void setPhysicsData();

  /// assign start pointers for arrays 2D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

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

#endif // ShockFitting_CopyRoeValues2_hh

