// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BoundaryConnectivity_hh
#define ShockFitting_BoundaryConnectivity_hh

//--------------------------------------------------------------------------//

#include <cassert>
#include <cstddef>

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
    m_dim(0), m_nbBnd(0), m_bInfo(NULL), m_elementNode(NULL), m_elementPtr(NULL) {}
  
  /// constructor
  /// @param dim          spatial dimension (2 or 3)
  /// @param nbBoundaries number of individual boundary patches
  /// @param bInfo        array of size @see nbBoundaries x 2, storing consecutively 
  ///                     (1) the number of elements (faces) in each boundary 
  ///                     (2) a pointer to the first element in each boundary
  /// @param bNode        element-node connectivity
  /// @param bPtr         pointer to the start of the corresponding element  
  BoundaryConnectivity(const unsigned dim,
		       const unsigned nbBoundaries, 
		       int* bInfo,
		       int* bNode, 
		       int* bPtr) 
  {reset(dim, nbBoundaries, bInfo, bNode, bPtr);}
  
  /// destructor
  ~BoundaryConnectivity() {}
  
  /// reset the connectivity
  void reset(const unsigned dim, const unsigned nbBoundaries, int* bInfo, int* bNode, int* bPtr)
  {m_dim = dim; m_nbBnd = nbBoundaries; m_bInfo = bInfo; m_elementNode = bNode; m_elementPtr = bPtr;}
  
  /// get the number of elements
  unsigned getNbElems(unsigned iBnd) const {return m_bInfo[iBnd*2];}
  
  /// get the pointer to the start of the specified boundary
  unsigned getBndPtr(unsigned iBnd) const {return m_bInfo[iBnd*2+1];}
  
  /// get the connectivity table for the specified boundary
  /// @param iBnd  boundary ID
  int* getTable(unsigned iBnd) const {return &m_elementNode[getBndPtr(iBnd)*getElemStride()];}
  
  /// get the connectivity pointers for the specified boundary
  /// @param iBnd  boundary ID
  int* getPtrs(unsigned iBnd) const {return &m_elementPtr[getBndPtr(iBnd)];}
  
  /// get node ID in element
  int getNodeID(unsigned iBnd, unsigned iElem, unsigned iNode) const 
  {assert(iElem < getNbElems(iBnd)); return getTable(iBnd)[getPtrs(iBnd)[iElem] + iNode];}
  
  /// get number of nodes in element
  int getNbNodes(unsigned iBnd, unsigned iElem) const 
  {assert(iElem < getNbElems(iBnd)); return m_elementPtr[getBndPtr(iBnd+1) - getBndPtr(iBnd)];}
  
private: // helper functions
  
  /// get the number of nodes per element (face)
  int getElemStride() const {return 2*m_dim-2;}
  
private: // data
  
  /// spatial dimension
  unsigned m_dim;
  
  /// number of boundary patches
  unsigned m_nbBnd;
  
  /// array of size m_nbBnd x 2, storing consecutively 
  /// (1) the number of elements (faces) in each boundary 
  /// (2) a pointer to the first element in each boundary
  int* m_bInfo;
  
  /// element-DOF connectivity, i.e. a list of all DOF IDs for each element 
  int* m_elementNode;
  
  /// array storing pointers to the beginning of each element
  int* m_elementPtr;
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_BoundaryConnectivity_hh
