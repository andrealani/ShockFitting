// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CopyRoeValues1_hh
#define ShockFitting_CopyRoeValues1_hh

//--------------------------------------------------------------------------//

#include "Framework/CopyMaker.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CopyRoeValues1, whose task is to overwrite the 
/// zroe values belonging to the background grid with the zroe values belonging
/// the the shocked one.

class CopyRoeValues1 : public CopyMaker {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CopyRoeValues1(const std::string& objectName);

  /// Destructor
  virtual ~CopyRoeValues1();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// copy new zroe values in the old ones
  virtual void copy();

private: // helper functions

  /// assign variables used in CopyRoevalues to MeshData pattern
  void setMeshData();

  /// assign variables used in CopyRoevalues to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers for arrays 2D
  void setAddress();

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
  std::vector <int>* M12M0;

  /// dummy variable to store starting pointer of Array2D
  unsigned start;
};

//--------------------------------------------------------------------------//

} // namespace Shockfitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_CopyRoeValues1_hh
