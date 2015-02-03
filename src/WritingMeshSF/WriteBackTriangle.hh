// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteBackTriangle_hh
#define ShockFitting_WriteBackTriangle_hh

//--------------------------------------------------------------------------//

#include <stdio.h>
#include <vector>
#include "Framework/WritingMesh.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a WriteBackTriangle, whose task is to write the
/// updates value of ZROE on triangle format file .node
/// File node format:
/// First line: <# of vertices> <dimension (must be 2)> <# of
/// attributes> <# of boundary markers (0 or 1)>
/// Remaining lines: <vertex #> <x> <y> [attributes]
/// [boundary marker] 

class WriteBackTriangle : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteBackTriangle(const std::string& objectName);

  /// Destructor
  virtual ~WriteBackTriangle();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// Write Triangle fomat file
  virtual void write();

private: // helper functions

  /// assign variables used in WriteBackTriangle to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteBackTriangle to PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointers for Array2D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shock points
  std::vector<unsigned>* npoin;

  /// mesh points state
  std::vector<double>* zroeVect;

  /// mesh points coordinates
  std::vector<double>* coorVect;

  /// code characterizing mesh points
  std::vector<int>* nodcod;

  /// name of the current output file
  std::string* fnameBack;

  /// mesh points state (in array storing)
  Array2D<double>* Zroe;

  /// mesh points coordinates (in array storing)
  Array2D<double>* XY;

  /// variable writing on node file
  FILE* file;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteBackTriangle
