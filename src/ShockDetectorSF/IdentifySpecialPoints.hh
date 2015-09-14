// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_IdentifySpecialPoints_hh
#define ShockFitting_IdentifySpecialPoints_hh

//--------------------------------------------------------------------------//

#include <sstream>

#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a IdentifySpecialPoints, whose task is to identify
/// the shock special points within the domain

/// @author Valentina De Amicis

class IdentifySpecialPoints {
public:

  /// Constructor
  IdentifySpecialPoints();

  /// Destructor
  ~IdentifySpecialPoints();

  /// set the special points
  void identify();  

private: // helper functions

  /// return the name of the class
  std::string getClassName() const { return "IdentifySpecialPoints"; }

  /// resize vectors and array
  void setSize();

  /// assign starting pointers for the array
  void setAddress();

  /// assign variables used in IdentifySpecialPoints to MeshData
  void setMeshData();

  /// assign variables used in IdentifySpecialPoints to PhysicsData
  void setPhysicsData();

  /// de-allocate the dynamic array
  void freeArray();

private: // data

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// boundary conditions string map
  std::vector<std::string>* BCmap;

  /// number of shock points foe each shock
  std::vector<unsigned>* nShockPoints;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// mesh boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// mesh points coordinates (in array storing)
  Array2D<double>* XY;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// file storing information
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_IdentifySpecialPoints
