// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BndryFacePtrFreez_hh
#define ShockFitting_BndryFacePtrFreez_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/Remeshing.hh"
#include "MathTools/Array2D.hh"
#include "Framework/FileLogManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a BndryFacePtrFreez, whose task is to set boundary faces
/// using the edge pointer. Skip interior constrained faces flag=999
/// This class is used when the starting grid is structured

class BndryFacePtrFreez: public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  BndryFacePtrFreez(const std::string& objectName);

  /// Destructor
  virtual ~BndryFacePtrFreez();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const { return "BndryFacePtrFreez"; }

  /// assign variables used in BdryFacePtr to MeshData pattern
  void setMeshData();

  /// assign array starting pointers
  void setAddress();

private: // data

  /// number of vertex
  unsigned* nvt;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// number of edges
  std::vector<unsigned>* nedge;

  /// vector characterizing boundary segments (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing edges (assignable to MeshData)
  std::vector<int>* edgptrVect;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// edgptr(0)(i-edge) 1° endpoint of i-edge
  /// edgptr(1)(i-edge) 2° endpoint of i-edge
  /// edgptr(2)(i-edge) boundary marker of i-edge
  Array2D <int>* edgptr;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_BndryFacePtrFreez_hh
