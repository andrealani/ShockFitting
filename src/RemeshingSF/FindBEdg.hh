// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FindBEdg_hh
#define ShockFitting_FindBEdg_hh

//----------------------------------------------------------------------------//

#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines FindBEdg, whose task is find the boundary edge
/// (of the background mesh) the shokc point (xsh, ysh) belongs to

class FindBEdg {
public:

  /// Constructor
  FindBEdg();

  /// Destructor
  ~FindBEdg();

  /// return findbedg
  int getBEdg(double, double);

  /// return s value
  double getS() const {return s;}

private: // helper functions

  /// assign values used in FindBEdg to MeshData pattern
  void setMeshData();

  /// assign values used in FindBEdg to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

private: // data

  /// value computed by class
  double s;

  /// space dimension
  unsigned* ndim;

  /// max number of shocks
  unsigned* nshmax;
  
  /// max number of shock edges
  unsigned* neshmax;
  
  /// number of boundar faces
  std::vector<unsigned>* nbfac;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// mesh points coordinates
  std::vector <double>* coorVect;

  /// vector characterizing boundary faces
  std::vector<int>* bndfacVect;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// mesh points coordinates (in array storing) 
  Array2D <double>* XY;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_FindBEdg_hh
  
