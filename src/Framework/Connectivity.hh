// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Connectivity_hh
#define ShockFitting_Connectivity_hh

//--------------------------------------------------------------------------//

#include <cassert>
#include <cstddef>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a basic connectivity object
/// 
/// @author Andrea Lani

class Connectivity {
public:
  /// default constructor
  Connectivity() : m_nbElems(0), m_elementNode(NULL), m_elementPtr(NULL) {}
  
  /// constructor
  Connectivity(const unsigned nbElements, int* eNode, int* ePtr) 
  {reset(nbElements, eNode, ePtr);}
  
  /// destructor
  ~Connectivity() {}
  
  /// reset the connectivity
  void reset(const unsigned nbElements, int* eNode, int* ePtr)
  {m_nbElems = nbElements; m_elementNode = eNode; m_elementPtr = ePtr;}
  
  /// get the number of elements
  unsigned getNbElems() const {return m_nbElems;}
  
  /// get the connectivity table
  int* getTable() const {return m_elementNode;}
  
  /// get the connectivity pointers
  int* getPtrs() const {return m_elementPtr;}
  
  /// get node ID in element
  int getNodeID(unsigned iElem, unsigned iNode) const 
  {assert(iElem < m_nbElems); return m_elementNode[m_elementPtr[iElem] + iNode];}
  
  /// get number of nodes in element
  int getNbNodes(unsigned iElem) const 
  {assert(iElem < m_nbElems); return m_elementPtr[iElem+1] - m_elementPtr[iElem];}
  
private:
  
  /// number of elements in the connectivity
  unsigned m_nbElems;
  
  /// element-DOF connectivity, i.e. a list of all DOF IDs for each element 
  int* m_elementNode;
  
  /// array storing pointers to the beginning of each element
  int* m_elementPtr;
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_Connectivity_hh
