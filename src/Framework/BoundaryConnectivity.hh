// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BoundaryConnectivity_hh
#define ShockFitting_BoundaryConnectivity_hh

//--------------------------------------------------------------------------//

#include <iostream>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <stdio.h>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a boundary connectivity object
/// 
/// @author Andrea Lani

class BoundaryConnectivity {
public:
  /// default constructor
  BoundaryConnectivity() : 
    m_dim(0), m_nbBnd(0), m_bInfo(NULL), m_faceNode(NULL), m_facePtr(NULL) {}
  
  /// constructor
  /// @param dim          spatial dimension (2 or 3)
  /// @param nbBoundaries number of individual boundary patches
  /// @param bInfo        array of size @see nbBoundaries x 2, storing consecutively 
  ///                     (1) the number of elements (faces) in each boundary 
  ///                     (2) a pointer to the first element in each boundary
  /// @param bNode        element-node connectivity
  /// @param bPtr         pointer to the start of the corresponding element  
  /// @param bName        names of the boundary patches
  BoundaryConnectivity(const unsigned dim,
		       const unsigned nbBoundaries, 
		       int* bInfo,
		       int* bNode, 
		       int* bPtr,
                       std::string* bName) 
  {reset(dim, nbBoundaries, bInfo, bNode, bPtr, bName);}
  
  /// destructor
  ~BoundaryConnectivity() {}
  
  /// reset the connectivity
  void reset(const unsigned dim, const unsigned nbBoundaries, int* bInfo, int* bNode, int* bPtr, std::string* bName)
  {m_dim = dim; m_nbBnd = nbBoundaries; m_bInfo = bInfo; m_faceNode = bNode; m_facePtr = bPtr; m_bName = bName;}

  /// get the number of boundary patches
  unsigned getNbBoundaries() const {return m_nbBnd;}
  
  /// get the number of faces
  unsigned getNbFaces(unsigned iBnd) const {return m_bInfo[iBnd*2];}
  
  /// get the pointer to the start of the specified boundary
  unsigned getBndPtr(unsigned iBnd) const {return m_bInfo[iBnd*2+1];}
  
  /// get the connectivity table for the specified boundary
  /// @param iBnd  boundary ID
  int* getTable(unsigned iBnd) const {return &m_faceNode[getBndPtr(iBnd)*getElemStride()];}

  /// get the connectivity pointers for the specified boundary
  /// @param iBnd  boundary ID
  int* getPtrs(unsigned iBnd) const {return &m_facePtr[getBndPtr(iBnd)];}
  
  /// get node ID in face 
  int getNodeID(unsigned iBnd, unsigned iFace, unsigned iNode) const 
  {assert(iFace < getNbFaces(iBnd)); return getTable(iBnd)[getPtrs(iBnd)[iFace] + iNode];}
  
  /// get number of nodes in Faces 
  int getNbNodes(unsigned iBnd, unsigned iFace) const 
  {assert(iFace < getNbFaces(iBnd)); return m_facePtr[getBndPtr(iBnd+1) - getBndPtr(iBnd)];}

  /// get name of the boundary patch
  /// @param iBnd  boundary ID
  std::string getBndName(unsigned iBnd) 
  {return m_bName[iBnd];}

private: // helper functions
  
  /// get the number of nodes per face
  int getElemStride() const {return 2*m_dim-2;}
  
private: // data
  
  /// spatial dimension
  unsigned m_dim;
  
  /// number of boundary patches
  unsigned m_nbBnd;
  
  /// array of size m_nbBnd x 2, storing consecutively 
  /// (1) the number of faces in each boundary 
  /// (2) a pointer to the first element in each boundary
  int* m_bInfo;
  
  /// face-DOF connectivity, i.e. a list of all DOF IDs for each face
  int* m_faceNode;
  
  /// array storing pointers to the beginning of each face 
  int* m_facePtr;

  /// array storing boundary patches names
  std::string* m_bName;
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_BoundaryConnectivity_hh
