// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ResetIniSol_hh
#define ShockFitting_ResetIniSol_hh

//--------------------------------------------------------------------------//

#include "Framework/WritingMesh.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ResetIniSol, whose task is 

class ResetIniSol : public WritingMesh {
public:

  /// Constructor
  ResetIniSol();

  /// Destructor
  virtual ~ResetIniSol();

  /// Set up object before its first use
  virtual void setup();

  /// Unset up object before its last use
  virtual void unsetup();

  /// reset initial solution
  virtual void write();

private: // helper function

  /// assign variables used in ResetIniSol to MeshData pattern
  void setMeshData();

  /// assign variables used in ResetIniSol to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

private: // data

  /// space dimension
  unsigned* ndim;

  /// number of degrees of freedom
  unsigned* ndof;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// number of mesh points
  unsigned* npoin;

  /// mesh status
  std::vector <double>* zroe;

  /// mesh points coordinates
  std::vector<double>* coor;

  /// mesh points coordinates (in array storing)
  Array2D <double>* w_XY;

  /// upstream shock points coordinates
  Array3D <double>* w_XYShu;

  /// downstream shock points coordinates
  Array3D <double>* w_XYShd;

  /// upstream shock points status
  Array3D <double>* w_ZRoeShu;

  /// downstream shock points status
  Array3D <double>* w_ZRoeShd;
};

//--------------------------------------------------------------------------//

}  // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ResetIniSol_hh
