// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReadCFmesh_hh
#define ShockFitting_ReadCFmesh_hh

//--------------------------------------------------------------------------//

#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/FileProcessing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ReadCFmesh, whose task is to read a CFmesh file 
/// (i.e. COOLFluiD's mesh+solution format) and store its data.
/// 
/// @author Andrea Lani
  
class ReadCFmesh : public FileProcessing {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  ReadCFmesh(const std::string& objectName);
  
  /// Destructor
  virtual ~ReadCFmesh();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the file processing 
  virtual void process();
    
  /// Get the solution field
  virtual void getSolutionField(Field* field)
  {
    Connectivity conn(m_totNbElem, &m_eState[0], &m_eptrs[0]); 
    field->reset(m_totNbStates, m_dimension + m_nbEqs + m_sizeExtraVars, conn, &m_states[0]);
  }

  /// Get the mesh field
  virtual void getMeshField(Field* field)
  {
    Connectivity conn(m_totNbElem, &m_eNode[0], &m_eptrn[0]); 
    field->reset(m_totNbNodes, m_dimension, conn, &m_nodes[0]);
  }  
  
private: // helper functions
  
  /// Get the number of times the file has been read on end (one after the other)
  unsigned getReadCount() const {return m_readCount;}
  
  /// Definition of ElementTypeData class
  class ElementTypeData {
  public:
    std::string shape;
    unsigned nbElems;
    unsigned nbNodes;
    unsigned nbStates;
  };
  
  /// Sets m_mapString2Reader, that maps a given string to a corresponding
  /// reader function
  void setMapString2Readers();
 
  /// Read an entry in the .CFmesh file
  bool readString(std::ifstream& file);	
  
  /// Reads the space dimension
  void readCFVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readSvnVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readCFmeshVersion(std::ifstream& fin);

  /// Reads the space dimension
  void readDimension(std::ifstream& fin);

  /// Reads the number of equations
  void readNbEquations(std::ifstream& fin);

  /// Reads the number of nodes
  void readNbNodes(std::ifstream& fin);

  /// Reads the nb of dofs state tensors and initialize with them the dofs
  void readNbStates(std::ifstream& fin);

  /// Reads the nb of extra variables associated with states
  void readNbExtraStateVars(std::ifstream& fin);
  
  /// Reads the nb of extra variables
  void readNbExtraVars(std::ifstream& fin);
  
  /// Reads the names of the extra variables associated with states
  void readExtraStateVarNames(std::ifstream& fin);

  /// Reads the names of the extra variables
  void readExtraVarNames(std::ifstream& fin);

  /// Reads the strides of the extra variables associated with states
  void readExtraStateVarStrides(std::ifstream& fin);

  /// Reads the strides of the extra variables
  void readExtraVarStrides(std::ifstream& fin);
  
  /// Reads the extra variables associated with states
  void readExtraVars(std::ifstream& fin);

  /// Reads the nb of elements
  void readNbElements(std::ifstream& fin);

  /// Reads the nb of element types
  void readNbElementTypes(std::ifstream& fin);

  /// Reads the order of the polynomial representation of the geometry
  void readGeometricPolyOrder(std::ifstream& fin);
  
  /// Reads the order of the polynomial representation of the solution
  void readSolutionPolyOrder(std::ifstream& fin);
  
  /// Reads the Type of the polynomial representation of the geometry
  void readGeometricPolyType(std::ifstream& fin);
  
  /// Reads the Type of the polynomial representation of the solution
  void readSolutionPolyType(std::ifstream& fin);

  /// Reads the element types (CFGeoShape::Type)
  void readElementTypes(std::ifstream& fin);

  /// Reads the nb of elements per type
  void readNbElementsPerType(std::ifstream& fin);

  /// Reads the nb of nodes per type
  void readNbNodesPerType(std::ifstream& fin);

  /// Reads the nb of dofs per type
  void readNbStatesPerType(std::ifstream& fin);

  /// Reads the list of nodes
  void readNodeList(std::ifstream& fin);

  /// Reads the list of state tensors and initialize the dofs
  void readStateList(std::ifstream& fin);

  /// Reads the data concerning the elements
  void readElementList(std::ifstream& fin);

  /// Reads the number of topological region sets and initialize the vector
  /// that will contain the all the topological region sets
  /// @pre the element list has been already read
  /// @pre Connection has been already been and set
  /// @pre some topological region sets have been already constructed (INNER_CELLS, INNER_FACES)
  void readNbTRSs(std::ifstream& fin);

  /// Reads the name of the current topological region sets
  void readTRSName(std::ifstream& fin);

  /// Reads the number of topological regions in the current
  /// topological region  set
  void readNbTRs(std::ifstream& fin);

  /// Reads the number of geometric entities in each topological
  /// region of the current topological region set
  void readNbGeomEnts(std::ifstream& fin);

  /// Reads the type of geometric entity in the current topological
  /// region set
  /// @pre  the read string must be "Face", "Cell" (or "Edge" in the future)
  /// @post the read string is converted in the corresponding CFGeoEnt::Type
  ///       by the method m_getCFGeoEnt::Type()
  void readGeomType(std::ifstream& fin);

  /// Reads all the lists of geometric entities, using them to build the
  /// corresponding topological region.
  /// Once that all the topological regions belonging to the current topological
  /// region set have been built, the topological region set itself is built.
  void readGeomEntList(std::ifstream& fin);
    
  /// Check if the given shape is valid
  bool isValidShape(const std::string& shape) const
  {
    for (unsigned i = 0; i < m_cfshapes.size(); ++i) {if (m_cfshapes[i] == shape) return true;}
    return false;
  }
  
  /// get the input file name
  std::string getInputFile() const 
  {
    using namespace std;
    assert(m_inputFiles.size() == 1); 
    string name = m_inputFiles[0];
    unsigned start = name.find(".CFmesh");
    name.erase(start);
    return string(name + "-iter_" + SConfig::to_str(m_iter) + ".CFmesh");
  }
  
private: // data
  
  /// Definition of geometric entity class
  struct GE {
    std::vector<unsigned> nodes;
    std::vector<unsigned> states;
  };
  
  /// Definition of topological region class
  struct TR {
    std::vector<GE> geo;
  };
  
  /// Definition of topological region set class
  struct TRS {
    std::string name;
    std::vector<TR> tr;
  };
  
  /// pointer to ReaderFunction
  typedef void (ReadCFmesh::*ReaderFunction)(std::ifstream& fin);
  
  /// type that maps a string read in the File with a ReaderFunction
  typedef std::map<std::string, 
		   ReaderFunction,
                   std::less<std::string> > MapString2Reader;
  
  /// map each string with a corresponding pointer to member function
  MapString2Reader m_mapString2Reader;
  
  /// space dimension
  unsigned m_dimension;
  
  /// number of equations
  unsigned m_nbEqs;
  
  /// total number of nodes
  unsigned m_totNbNodes;  
  
  /// total number of states
  unsigned m_totNbStates;
  
  /// total number of elements
  unsigned m_totNbElem;
  
  /// total number of element types
  unsigned m_totNbElemTypes;
  
  /// number of extra variables (redundant)
  unsigned m_nbExtraVars;
  
  /// number of extra state variables
  unsigned m_nbExtraStateVars;
  
  /// total number of extra variables
  unsigned m_sizeExtraVars;
  
  /// list of supported CFmesh shapes
  std::vector<std::string> m_cfshapes;
  
  /// list of extra variable names (redundant)
  std::vector<std::string> m_extraVarNames;
  
  /// list of extra state variable names
  std::vector<std::string> m_extraStateVarNames;
  
  /// list of extra variable strides (redundant)
  std::vector<unsigned> m_extraVarStrides;
  
  /// list of extra state variable strids
  std::vector<unsigned> m_extraStateVarStrides;
  
  /// list of element types data 
  std::vector<ElementTypeData> m_elemTypes;
  
  /// storage of coordinates
  std::vector<double> m_nodes;
  
  /// storage of state values
  std::vector<double> m_states;
  
  /// element-node connectivity
  std::vector<int> m_eNode;
  
  /// element-state connectivity
  std::vector<int> m_eState;
  
  /// element-node pointers
  std::vector<int> m_eptrn;
  
  /// element-state pointers
  std::vector<int> m_eptrs;
  
  /// list of TRSs
  std::vector<TRS> m_trs;
  
  /// current ID of TRS
  int m_currTrsID;
  
  /// current ID of TR
  int m_currTrID;
  
  /// reader count
  unsigned m_readCount;
  
  /// current iteration number
  unsigned m_iter;
  
  /// location at which the state list starts
  long  m_startStateList;
  
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_ReadCFmesh_hh
