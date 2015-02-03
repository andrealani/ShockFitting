// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ChangeBndryPtr_hh
#define ShockFitting_ChangeBndryPtr_hh

//----------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a ChangeBndryPtr, whose task is to update the 
/// boundary data structure to take into account the boundary phantom nodes

class ChangeBndryPtr : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ChangeBndryPtr(const std::string& objectName);

  /// Destructor
  virtual ~ChangeBndryPtr();

  /// Setup object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// change boundary node pointer
  virtual void remesh();

private: // helper functions

  /// search the i-node IPOIN in inodptr vector if the node has been
  /// dis-actived (nodcod=-2)
  void lookForNode(int);

  /// remove the i-node
  void removeNode(int);

  /// update the i-face=nodptr(ipos,2)
  void updateIface();

  /// remove the i-face=nodptr(ipos,3)
  void removeIface();

  /// return class name
  std::string getClassName() {return std::string("ChangeBndryPtr");}

  /// assign starting pointers to arrays
  void setAddress();

  /// assign values used in ChangeBndryPtr to MeshData pattern
  void setMeshData();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of points in the mesh
  std::vector<unsigned>* npoin;

  /// number of boundary points
  std::vector<unsigned>* nbpoin;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// vector characterizing boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing boundary points (assignable to MeshData)
  std::vector<int>* nodptrVect;

  /// nodptr(i-poin)(0) global code number
  /// nodptr(i-poin)(1) one of the edges it belongs to
  /// nodptr(i-poin)(2) other edge
  Array2D <int>* nodptr;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// position of i-node in nodptr vector
  int ipos;

  /// i-face
  int iface;

  /// dummy vector used to store bndfac elements
  std::vector <int> inode;

  /// store file log info
  FileLogManip logfile;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ChangeBndryPtr_hh
