// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <limits>
#include <valarray>

#include "FileProcessingSF/ReadCFmesh.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Clist.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReadCFmesh, FileProcessing> readCFmeshProv("ReadCFmesh"); 
  
//--------------------------------------------------------------------------//
  
ReadCFmesh::ReadCFmesh(const std::string& objectName) :
  FileProcessing(objectName), 
  m_mapString2Reader(),
  m_dimension(0),
  m_nbEqs(0),
  m_totNbNodes(0),
  m_totNbStates(0),
  m_totNbElem(0),
  m_totNbElemTypes(0),
  m_nbExtraVars(0),
  m_nbExtraStateVars(0),
  m_sizeExtraVars(),
  m_cfshapes(),
  m_extraVarNames(),
  m_extraStateVarNames(),
  m_extraVarStrides(),
  m_extraStateVarStrides(),
  m_elemTypes(),
  m_nodes(),
  m_states(),
  m_eNode(),
  m_eState(),
  m_eptrn(),
  m_eptrs(),
  m_trs(),
  m_currTrsID(-1),
  m_currTrID(-1),
  m_readCount(0),
  m_iter(0),
  m_startStateList(0)
{
}
 
//--------------------------------------------------------------------------//

ReadCFmesh::~ReadCFmesh()
{
}

//--------------------------------------------------------------------------//

void ReadCFmesh::setup()
{  
  LogToScreen(VERBOSE, "ReadCFmesh::setup() => start\n");

  setMapString2Readers();
  
  /// store CFmesh element shapes
  m_cfshapes.reserve(7);
  m_cfshapes.push_back("Line");
  m_cfshapes.push_back("Triag");
  m_cfshapes.push_back("Quad");
  m_cfshapes.push_back("Tetra");
  m_cfshapes.push_back("Pyram");
  m_cfshapes.push_back("Prism");
  m_cfshapes.push_back("Hexa");
  
  // initialize the local iteration
  m_iter = getStartIter(); 
  
  LogToScreen(VERBOSE, "ReadCFmesh::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::unsetup() 
{ 
  LogToScreen(VERBOSE, "ReadCFmesh::unsetup()\n");
} 
 
//--------------------------------------------------------------------------// 

void ReadCFmesh::process() 
{  
  LogToScreen(INFO, "ReadCFmesh::process()\n");
  
  const long MAXLINESINFILE = std::numeric_limits<long>::max()-1;
    
  if (!fileExists(getInputFile().c_str())) {
    LogToScreen(INFO,  "ReadCFmesh::process() => ERROR: file <" << getInputFile() << "> does not exist in current directory!\n"); abort(); 
  }
  
  ifstream file(getInputFile().c_str());
  bool keepOnReading = true;
  long linesRead = 0;
  
  do {
    keepOnReading = readString(file);
    
    if (++linesRead > MAXLINESINFILE) {
      LogToScreen(INFO,  "ReadCFmesh::process() => File too long, probably misformatted."
		  "\nSee void FileReader::readFromFile for more info\n");
    }
  } while (keepOnReading);
  
  file.close();
  
  m_readCount++;
  
  m_iter += getIORate();
  
  /// rescale the field variables 
  Field f; 
  getSolutionField(&f);
  rescaleField(&f);
}

//--------------------------------------------------------------------------// 

void ReadCFmesh::setMapString2Readers()
{
  m_mapString2Reader["!COOLFLUID_VERSION"]     = &ReadCFmesh::readCFVersion;
  m_mapString2Reader["!COOLFLUID_SVNVERSION"]  = &ReadCFmesh::readSvnVersion;
  m_mapString2Reader["!CFMESH_FORMAT_VERSION"] = &ReadCFmesh::readCFmeshVersion;
  m_mapString2Reader["!NB_DIM"]             = &ReadCFmesh::readDimension;
  m_mapString2Reader["!NB_EQ"]              = &ReadCFmesh::readNbEquations;
  m_mapString2Reader["!NB_NODES"]           = &ReadCFmesh::readNbNodes;
  m_mapString2Reader["!NB_STATES"]          = &ReadCFmesh::readNbStates;
  m_mapString2Reader["!NB_EXTRA_SVARS"]     = &ReadCFmesh::readNbExtraStateVars;
  m_mapString2Reader["!NB_EXTRA_VARS"]     = &ReadCFmesh::readNbExtraVars;
  m_mapString2Reader["!EXTRA_SVARS_NAMES"]  = &ReadCFmesh::readExtraStateVarNames;
  m_mapString2Reader["!EXTRA_VARS_NAMES"]  = &ReadCFmesh::readExtraVarNames;
  m_mapString2Reader["!EXTRA_SVARS_STRIDES"]  = &ReadCFmesh::readExtraStateVarStrides;
  m_mapString2Reader["!EXTRA_VARS_STRIDES"]  = &ReadCFmesh::readExtraVarStrides;
  m_mapString2Reader["!EXTRA_VARS"]         = &ReadCFmesh::readExtraVars;
  m_mapString2Reader["!NB_ELEM"]            = &ReadCFmesh::readNbElements;
  m_mapString2Reader["!NB_ELEM_TYPES"]      = &ReadCFmesh::readNbElementTypes;
  m_mapString2Reader["!GEOM_POLYTYPE"]      = &ReadCFmesh::readGeometricPolyType;
  m_mapString2Reader["!SOL_POLYTYPE"]       = &ReadCFmesh::readSolutionPolyType;
  m_mapString2Reader["!GEOM_POLYORDER"]     = &ReadCFmesh::readGeometricPolyOrder;
  m_mapString2Reader["!SOL_POLYORDER"]      = &ReadCFmesh::readSolutionPolyOrder;
  m_mapString2Reader["!LIST_NODE"]          = &ReadCFmesh::readNodeList;
  m_mapString2Reader["!LIST_STATE"]         = &ReadCFmesh::readStateList;
  m_mapString2Reader["!NB_TRSs"]            = &ReadCFmesh::readNbTRSs;
  m_mapString2Reader["!TRS_NAME"]           = &ReadCFmesh::readTRSName;
  m_mapString2Reader["!NB_TRs"]             = &ReadCFmesh::readNbTRs;
  m_mapString2Reader["!NB_GEOM_ENTS"]       = &ReadCFmesh::readNbGeomEnts;
  m_mapString2Reader["!GEOM_TYPE"]          = &ReadCFmesh::readGeomType;
  m_mapString2Reader["!LIST_GEOM_ENT"]      = &ReadCFmesh::readGeomEntList;
  m_mapString2Reader["!ELEM_TYPES"]         = &ReadCFmesh::readElementTypes;
  m_mapString2Reader["!NB_ELEM_PER_TYPE"]   = &ReadCFmesh::readNbElementsPerType;
  m_mapString2Reader["!NB_NODES_PER_TYPE"]  = &ReadCFmesh::readNbNodesPerType;
  m_mapString2Reader["!NB_STATES_PER_TYPE"] = &ReadCFmesh::readNbStatesPerType;
  m_mapString2Reader["!LIST_ELEM"]          = &ReadCFmesh::readElementList;
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readCFVersion(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readCFVersion() start\n");
  std::string version; fin >> version;
  LogToScreen(DEBUG_MIN,   "CF version : " << version << "\n");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readCFVersion() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readSvnVersion(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSvnVersion() start\n");
  std::string version; fin >> version;
  LogToScreen(DEBUG_MIN,   "svn version : " << version << "\n");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSvnVersion() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readCFmeshVersion(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readCFmeshVersion() start\n");
  std::string version; fin >> version;
  LogToScreen(DEBUG_MIN,   "CFmesh version : " << version << "\n");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readCFmeshVersion() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readDimension(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readDimension() start\n");
  fin >> m_dimension;
  validate ((m_dimension == 2 || m_dimension == 3),
	    string("Invalid dimension (" + to_str(m_dimension) + ") in CFmesh"));
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readDimension() end\n");
}
  
//--------------------------------------------------------------------------//
  
void ReadCFmesh::readNbEquations(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbEquations() start\n");
  fin >> m_nbEqs;
  
  validate (m_nbEqs > 0, 
	    string("Invalid nb of equations (" + to_str(m_dimension) + ") in CFmesh"));
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbEquations() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbNodes(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbNodes() start\n");
  
  int nbNonUpdatableNodes = 0;
  fin >> m_totNbNodes >> nbNonUpdatableNodes;
  
  validate (m_totNbNodes > 0, 
	    string("Invalid nb of nodes (" + to_str(m_totNbNodes) + ") in CFmesh"));
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbNodes() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbStates(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbStates() start\n");
  int nbNonUpdatableStates = 0;
  fin >> m_totNbStates >> nbNonUpdatableStates;
  
  validate (m_totNbStates > 0, 
	    string("Invalid nb of states (" + to_str(m_totNbStates) + ") in CFmesh"));
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbStates() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbElements(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbElements() start\n");
  fin >> m_totNbElem;
  
  validate (m_totNbElem > 0,
	    string("Invalid nb of elements (" + to_str(m_totNbElem) + ") in CFmesh"));
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbElements() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbElementTypes(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbElementTypes() start\n");
  fin >> m_totNbElemTypes;
  
  validate (m_totNbElemTypes > 0 && m_totNbElemTypes < 5, 
	    string("Invalid nb of element types (" + to_str(m_totNbElemTypes) + ") in CFmesh"));
  
  m_elemTypes.resize(m_totNbElemTypes);
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbElementTypes() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readGeometricPolyOrder(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeometricPolyOrder() start\n");
  unsigned order = 0; fin >> order;
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeometricPolyOrder() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readSolutionPolyOrder(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSolutionPolyOrder() start\n");
  unsigned order = 0; fin >> order;
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSolutionPolyOrder() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readGeometricPolyType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeometricPolyType() start\n");
  unsigned pType = 0; fin >> pType;
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeometricPolyType() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readSolutionPolyType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSolutionPolyType() start\n");
  unsigned sType = 0; fin >> sType;
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readSolutionPolyType() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readElementTypes(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementTypes() start\n");
  assert(m_totNbElemTypes > 0);
  for (unsigned i = 0; i < m_totNbElemTypes; ++i) {
    fin >> m_elemTypes[i].shape;
    
    validate (isValidShape(m_elemTypes[i].shape),
	      string("Invalid shape (" + m_elemTypes[i].shape + ") in CFmesh"));
  }
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementTypes() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbElementsPerType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementTypes() start\n");
  unsigned sumNbElems = 0;
  for (unsigned i = 0; i < m_totNbElemTypes; ++i) {
    fin >> m_elemTypes[i].nbElems;
    sumNbElems += m_elemTypes[i].nbElems;
  }
  
  validate (sumNbElems == m_totNbElem, 
	    string("Invalid sum of elements (" + to_str(sumNbElems) + " != " + to_str(m_totNbElem) + ") in CFmesh"));
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementTypes() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbNodesPerType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbNodesPerType() start\n");
  unsigned nbNPT = 0;
  for (unsigned i = 0; i < m_totNbElemTypes; ++i) {
    fin >> nbNPT;
    m_elemTypes[i].nbNodes = nbNPT;
    
    /// @todo try to avoid having these hardcoded (only working for P1 geometry elements)
    validate ((m_dimension == 1 && nbNPT == 2) ||
	      (m_dimension == 2 && nbNPT == Clist<unsigned>(3,4,6,8,9,10)) ||
	      (m_dimension == 3 && nbNPT == Clist<unsigned>(4,5,6,8,10,13,14,15,18,20,27)),
	      string("Invalid nb of nodes (" + to_str(nbNPT) + ") per element " + 
		     m_elemTypes[i].shape + ") in CFmesh"));
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbNodesPerType() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbStatesPerType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbStatesPerType() start\n");
  
  unsigned nbSPT = 0;
  for (unsigned i = 0; i < m_totNbElemTypes; ++i) {
    fin >> nbSPT;
    m_elemTypes[i].nbStates = nbSPT;
    
    validate (nbSPT > 0,
	      string("Invalid nb of states (" + to_str(nbSPT) + ") per element " + 
		     m_elemTypes[i].shape + ") in CFmesh"));
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbStatesPerType() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNodeList(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNodeList() start\n");
  m_nodes.reserve(m_dimension*m_totNbNodes);
  
  double entry = 0.;
  for (unsigned i = 0; i < m_totNbNodes; ++i) {
    for (unsigned d = 0; d < m_dimension; ++d) {
      fin >> entry;
      m_nodes.push_back(entry);
    }
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNodeList() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readStateList(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readStateList() start\n");
  bool isWithSolution = false;
  fin >> isWithSolution;
  validate (isWithSolution, "Solution not present in CFmesh file!");
  
  m_sizeExtraVars = 0;
  if (m_nbExtraStateVars > 0) {
    m_sizeExtraVars = std::accumulate(m_extraStateVarStrides.begin(), m_extraStateVarStrides.end(), 0);
  }
  m_states.reserve((m_dimension + m_nbEqs + m_sizeExtraVars)*m_totNbStates);
  
  Connectivity elements(m_totNbElem, &m_eNode[0], &m_eptrn[0]);
  valarray<double> centroid(m_dimension);
  double entry = 0.;
  
  for (unsigned i = 0; i < m_totNbStates; ++i) {
    // x,y,z coordinates are computed on the fly while reading the states 
    if (m_totNbNodes == m_totNbStates) {
      // vertex centered meshes
      const unsigned start = i*m_dimension;
      for (unsigned d = 0; d < m_dimension; ++d) {
	m_states.push_back(m_nodes[start + d]);
      }
    }
    else {
      // cell-centered meshes: position vector coincides with cell centroid
      centroid = 0.;
      const unsigned nbElemNodes = elements.getNbNodes(i);
      validate(nbElemNodes > 1 && nbElemNodes <= 8, "ReadCFmesh::readStateList() => wrong nbElemNodes!");
      for (unsigned iNode = 0; iNode < nbElemNodes; ++iNode) {
	const unsigned nodeID = elements.getNodeID(i, iNode);
	const unsigned start = nodeID*m_dimension;
	for (unsigned d = 0; d < m_dimension; ++d) {
	  centroid[d] += m_nodes[start+d];
	}
      }
      centroid /= nbElemNodes;
      for (unsigned d = 0; d < m_dimension; ++d) {
	m_states.push_back(centroid[d]);
      }
    }
    
    // state vector
    for (unsigned d = 0; d < m_nbEqs; ++d) {
      fin >> entry;
      m_states.push_back(entry);
    }
    // extra state variables
    for (unsigned d = 0; d < m_nbExtraStateVars; ++d) {
      const unsigned nbExs =  m_extraStateVarStrides[d];
      for (unsigned s = 0; s < nbExs; ++s) {
	fin >> entry;
	m_states.push_back(entry);
      }
    }
  }
    
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readStateList() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readElementList(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementList() start\n");
  
  unsigned elemID = 0;
  unsigned sizeElemNodeVec  = 0;
  unsigned sizeElemStateVec = 0;
  for (unsigned ie = 0; ie < m_elemTypes.size(); ++ie) {
    const ElementTypeData& el = m_elemTypes[ie];
    sizeElemNodeVec  += el.nbElems*el.nbNodes;
    sizeElemStateVec += el.nbElems*el.nbStates;
    elemID += el.nbElems;
  }
  validate(elemID == m_totNbElem, "ReadCFmesh::readElementList() => wrong nb elements");
  
  m_eNode.resize(sizeElemNodeVec);
  m_eState.resize(sizeElemStateVec);
  m_eptrn.resize(m_totNbElem + 1);
  m_eptrs.resize(m_totNbElem + 1);
  
  unsigned nbEEminProc       = m_totNbElem;
  unsigned nbEEminProcMinus1 = nbEEminProc - 1;
  unsigned ncount = 0;
  unsigned scount = 0;
  unsigned ipos = 0;
  unsigned iElemBegin = 0;
  for (unsigned iType = 0; iType < m_totNbElemTypes; ++iType) {
    const unsigned nbNodesInElem  = m_elemTypes[iType].nbNodes;
    const unsigned nbStatesInElem = m_elemTypes[iType].nbStates;
    const unsigned nbElementsPerType = m_elemTypes[iType].nbElems;
    const unsigned iElemEnd = iElemBegin + nbElementsPerType;
    
    // loop over the elements in this type
    for (unsigned iElem = iElemBegin; iElem < iElemEnd; ++iElem) {
      m_eptrn[ipos] = ncount;
      m_eptrs[ipos] = scount;
      
      for (unsigned j = 0; j < nbNodesInElem; ++j, ++ncount) {
	fin >> m_eNode[ncount];
	validate(m_eNode[ncount] < (int)m_totNbNodes, string("Invalid node ID (" + to_str(m_eNode[ncount]) + ")"));
      }
      for (unsigned j = 0; j < nbStatesInElem; ++j, ++scount) {
	fin >> m_eState[scount];
	validate(m_eState[scount] < (int)m_totNbStates, string("Invalid state ID (" + to_str(m_eState[scount]) + ")"));
      }
      
      if (ipos == nbEEminProcMinus1) {
	m_eptrn[nbEEminProc] = m_eptrn[ipos] + nbNodesInElem;
	m_eptrs[nbEEminProc] = m_eptrs[ipos] + nbStatesInElem;
      }
      ++ipos;
    }
    
    iElemBegin +=  nbElementsPerType;
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readElementList() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbTRSs(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbTRSs() start\n");

  unsigned nbTRSs = 0;
  fin >> nbTRSs;
  validate (nbTRSs > 0, "Number of TRSs in file must be at least 1");
  // WATCH OUT: we do not account for TRS reduction as in COOLFluiD's CFmesh reader
  
  // resize the storage of TRSs to be filled later
  m_trs.resize(nbTRSs);
  m_currTrsID = -1; // initialize counter for TRSs
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbTRSs() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readTRSName(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readTRSName() start\n");
  m_currTrsID++;   // increment counter for TRS (a new one starts)
  m_currTrID = -1; // initialize counter for TRs within the current TRS
  
  fin >> m_trs[m_currTrsID].name;
  LogToScreen(DEBUG_MIN,  "ReadCFmesh::readTRSName() => TRS[" << m_currTrsID << "].name is " 
       << m_trs[m_currTrsID].name  << "\n");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readTRSName() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbTRs(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbTRs() start\n");
  unsigned nbTRsInTRS = 0;
  fin >> nbTRsInTRS;
  validate(nbTRsInTRS > 0, "Invalid nb of TRs in CFmesh : must be > 0");
  
  m_trs[m_currTrsID].tr.resize(nbTRsInTRS);
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbTRs() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbGeomEnts(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbGeomEnts() start\n");
  const unsigned nbTRsInTRS = m_trs[m_currTrsID].tr.size();
  unsigned nbGEinTR = 0;
  for (unsigned i = 0; i < nbTRsInTRS; ++i) {
    fin >> nbGEinTR;
    validate(nbGEinTR > 0, "Invalid nb of GEs in TR in CFmesh : must be > 0");
    m_trs[m_currTrsID].tr[i].geo.resize(nbGEinTR);
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbGeomEnts() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readGeomType(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeomType() start\n");
  std::string geomName = "";
  fin >> geomName;
  validate(geomName == "Face", "GEOM_TYPE != Face not supported in CFmesh");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeomType() end\n");
}
	  
//--------------------------------------------------------------------------//

void ReadCFmesh::readGeomEntList(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeomEntList() start\n");
  const unsigned nbTRs = m_trs[m_currTrsID].tr.size();

  // loop only in the new TRs, which have not been read yet
  for (unsigned iTR = 0; iTR < nbTRs; ++iTR) {
    const unsigned nbTRGeos = m_trs[m_currTrsID].tr[iTR].geo.size();
    for (unsigned iGeo = 0; iGeo < nbTRGeos; ++iGeo) {
      unsigned nbNodesInGeo = 0;
      unsigned nbStatesInGeo = 0;
      fin >> nbNodesInGeo >> nbStatesInGeo;
      GE& geo = m_trs[m_currTrsID].tr[iTR].geo[iGeo];
      geo.nodes.resize(nbNodesInGeo);
      geo.states.resize(nbStatesInGeo);
      
      for(unsigned n = 0; n < nbNodesInGeo; ++n) {
        fin >> geo.nodes[n];
        assert(geo.nodes[n] < m_totNbNodes);
      }
      
      for(unsigned s = 0; s < nbStatesInGeo; ++s) {
        fin >> geo.states[s];
        assert(geo.states[s] < m_totNbStates);
      }
    } // loop iGeo
  }
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readGeomEntList() end\n");
}
	  
//--------------------------------------------------------------------------//

bool ReadCFmesh::readString(ifstream& file)
{
  std::string key = "";
  file >> key;
  
  LogToScreen(DEBUG_MIN,  "CFmesh key = " << key << "\n");
  
  // check end of file
  if (key != "!END") {

    MapString2Reader::iterator key_pair = m_mapString2Reader.find(key);

    // check if key exists else ignore it
    if (key_pair != m_mapString2Reader.end()) {
      
      if (getReadCount() == 0) {
	if (key == "!LIST_STATE") {m_startStateList = file.tellg();}
      }
      
      // if this is the second reading, only read list of states
      if (getReadCount() > 0) {
	// read state list
	file.seekg(m_startStateList);
	readStateList(file);
	return false;
      }
      else {
	LogToScreen(DEBUG_MIN,   "Read CFmesh Key: " << key << "\n");
	ReaderFunction function = key_pair->second;
	assert(function != NULL);
	(this->*function)(file);
      }
    }
    else {
      std::string msg = "Key in CFmesh is not valid:" + key;
      throw ConfigException(msg);
    }
    // keep reading
    return true;
  }
  
  // found end of file
  return false;
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readNbExtraVars(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbExtraVars() start\n");
  
  int nbExtraVars = 0;
  fin >> nbExtraVars;
  validate(nbExtraVars >= 0,  "Negative number of extra variables in CFmesh");
  m_nbExtraVars = nbExtraVars;
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbExtraVars() end\n");
}
      
      
//--------------------------------------------------------------------------//

void ReadCFmesh::readNbExtraStateVars(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbExtraStateVars() start\n");

  int nbExtraStateVars = 0;
  fin >> nbExtraStateVars;
  validate(nbExtraStateVars >= 0,  "Negative number of extra state variables in CFmesh");
  m_nbExtraStateVars = nbExtraStateVars;
  
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readNbExtraStateVars() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readExtraVarNames(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVarNames() start\n");
  m_extraVarNames.resize(m_nbExtraVars);
  for (unsigned i = 0; i < m_nbExtraVars; ++i) {
    fin >> m_extraVarNames[i];
  }
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVarNames() end\n");
}
      
//--------------------------------------------------------------------------//

void ReadCFmesh::readExtraStateVarNames(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraStateVarNames() start\n");
  m_extraStateVarNames.resize(m_nbExtraStateVars);
  for (unsigned i = 0; i < m_nbExtraStateVars; ++i) {
    fin >> m_extraStateVarNames[i];
  }
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraStateVarNames() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readExtraVarStrides(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVarStrides() start\n");
  m_extraVarStrides.resize(m_nbExtraVars);
  for (unsigned i = 0; i < m_nbExtraVars; ++i) {
    fin >> m_extraVarStrides[i];
  }
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVarStrides() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readExtraStateVarStrides(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraStateVarStrides() start\n");
  m_extraStateVarStrides.resize(m_nbExtraStateVars);
  for (unsigned i = 0; i < m_nbExtraStateVars; ++i) {
    fin >> m_extraStateVarStrides[i];
  }
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraStateVarStrides() end\n");
}

//--------------------------------------------------------------------------//

void ReadCFmesh::readExtraVars(ifstream& fin)
{
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVars() start\n");
  LogToScreen(DEBUG_MIN,   "ReadCFmesh::readExtraVars() end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
