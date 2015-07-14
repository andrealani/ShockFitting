// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>
#include "Framework/ShockFittingObj.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/BoundaryConnectivity.hh"
#include "Framework/Connectivity.hh"
#include "Framework/Field.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/ConfigFileReader.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ShockFittingObj, ShockFittingObj> 
ShockFittingObjProv("ShockFittingObj");

//--------------------------------------------------------------------------//
  
ShockFittingObj::ShockFittingObj(const std::string& objectName) : 
  Counter(),
  ConfigObject(objectName)
{
  m_fInterpolator = vector<PAIR_TYPE(FieldInterpolator)>();
  addOption("FieldInterpolatorList",&m_fInterpolator, 
	    "List of the names of field interpolators"); 
  
  m_fProcessing = vector<PAIR_TYPE(FileProcessing)>();
  addOption("FileProcessingList",&m_fProcessing, 
	    "List of the names of file processing"); 

  m_mGenerator = vector<PAIR_TYPE(MeshGenerator)>();
  addOption("MeshGeneratorList",&m_mGenerator,
            "List of the names of mesh generator");

  m_sDetector = vector<PAIR_TYPE(ShockDetector)>();
  addOption("ShockDetectorList",&m_sDetector,
            "List of the names of shock detectors");

  m_fRemeshing = vector<PAIR_TYPE(Remeshing)>();
  addOption("RemeshingList",&m_fRemeshing,
	    "List of the names of field remeshing");

  m_wMesh = vector<PAIR_TYPE(WritingMesh)>();
  addOption("WritingMeshList",&m_wMesh,
            "List of the names of objects writing");

  m_fConverter = vector<PAIR_TYPE(Converter)>();
  addOption("ConverterList",&m_fConverter,
            "List of the names of converter objects");

  m_cfdSolver = vector<PAIR_TYPE(CFDSolver)>();
  addOption("CFDSolverList",&m_cfdSolver,
            "List of the names of cfd solver objects");

  m_cMaker = vector<PAIR_TYPE(CopyMaker)>();
  addOption("CopyMakerList",&m_cMaker,
            "List of the names of copy maker objects");

  m_sUpdater = vector<PAIR_TYPE(StateUpdater)>();
  addOption("StateUpdaterList",&m_sUpdater,
            "List of the names of state updaters");
}
  
//--------------------------------------------------------------------------//

ShockFittingObj::~ShockFittingObj()
{
}

//--------------------------------------------------------------------------//

void ShockFittingObj::configure(SConfig::OptionMap& cmap,
				 const std::string& prefix)
{
  LogToScreen(VERBOSE, "ShockFittingObj::configure() => start\n");
  
  ConfigObject::configure(cmap, prefix);
  
  if (ConfigFileReader::isFirstConfig()) {
    
    // create the field interpolators
    createList<FieldInterpolator>(m_fInterpolator);
    
    // create the file processing
    createList<FileProcessing>(m_fProcessing);

    // create the mesh generator
    createList<MeshGenerator>(m_mGenerator);

    // create the shock detector
    createList<ShockDetector>(m_sDetector);

    // create the remeshing
    createList<Remeshing>(m_fRemeshing);

    // create the writing objects
    createList<WritingMesh>(m_wMesh);

    // create the converter objects
    createList<Converter>(m_fConverter);

    // create the object calling cfd solver 
    createList<CFDSolver>(m_cfdSolver);

    // create the copy maker objects
    createList<CopyMaker>(m_cMaker);

    // create the state updaters list
    createList<StateUpdater>(m_sUpdater);
  }

  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    configureDeps (cmap, m_fInterpolator[i].ptr().get());
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    configureDeps (cmap, m_fProcessing[i].ptr().get());
  }

  // configure physics data
  configureDeps (cmap, &PhysicsData::getInstance());
  
  // configure mesh data
  configureDeps (cmap, &MeshData::getInstance());

  // configure the mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    configureDeps (cmap, m_mGenerator[i].ptr().get());
  }

  // configure the shock detector
  for (unsigned i = 0; i < m_sDetector.size(); ++i) {
    configureDeps (cmap, m_sDetector[i].ptr().get());
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    configureDeps (cmap, m_fRemeshing[i].ptr().get());
  }

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
   configureDeps (cmap, m_wMesh[i].ptr().get());
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
   configureDeps (cmap, m_fConverter[i].ptr().get());
  }

  // configure the object calling CFD solver objects 
  for (unsigned i = 0; i < m_cfdSolver.size(); ++i) {
   configureDeps (cmap, m_cfdSolver[i].ptr().get());
  }

  // configure the copy maker objects
  for (unsigned i = 0; i < m_cMaker.size(); ++i) {
   configureDeps (cmap, m_cMaker[i].ptr().get());
  }

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
   configureDeps (cmap, m_sUpdater[i].ptr().get());
  }

  LogToScreen(VERBOSE, "ShockFittingObj::configure() => end\n");
}

//--------------------------------------------------------------------------//

void ShockFittingObj::setup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::setup() => start\n");
 
  createPhysicsData();

  createMeshData();

  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->setup();
  }

  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->setup();
  } 

  // configure the mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    m_mGenerator[i].ptr()->setup();
  }

  // configure the shock detector
  for (unsigned i = 0; i < m_sDetector.size(); ++i) {
    m_sDetector[i].ptr()->setup();
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    m_fRemeshing[i].ptr()->setup();
  }

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
    m_wMesh[i].ptr()->setup();
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
    m_fConverter[i].ptr()->setup();
  }

  // configure the object calling CFD solver objects
  for (unsigned i = 0; i < m_cfdSolver.size(); ++i) {
    m_cfdSolver[i].ptr()->setup();
  }

  // configure the copy maker objects
  for (unsigned i = 0; i < m_cMaker.size(); ++i) {
    m_cMaker[i].ptr()->setup();
  }

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
    m_sUpdater[i].ptr()->setup();
  }

  LogToScreen(VERBOSE, "ShockFittingObj::setup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::unsetup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => start\n");

  deleteMeshData();

  deletePhysicsData();
 
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->unsetup();
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->unsetup();
  }

  // configure the mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    m_mGenerator[i].ptr()->unsetup();
  }

  // configure the shock detector
  for (unsigned i = 0; i < m_sDetector.size(); ++i) {
    m_sDetector[i].ptr()->unsetup();
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    m_fRemeshing[i].ptr()->unsetup();
  }

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
    m_wMesh[i].ptr()->unsetup();
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
    m_fConverter[i].ptr()->unsetup();
  }

  // configure the object calling CFD solver objects
  for (unsigned i = 0; i < m_cfdSolver.size(); ++i) {
    m_cfdSolver[i].ptr()->unsetup();
  }

  // configure the copy maker objects
  for (unsigned i = 0; i < m_cMaker.size(); ++i) {
    m_cMaker[i].ptr()->unsetup();
  }

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
    m_sUpdater[i].ptr()->unsetup();
  }

  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::process()
{
  LogToScreen(VERBOSE, "ShockFittingObj::process()\n");
}

//--------------------------------------------------------------------------//
  
void ShockFittingObj::processField(const BoundaryConnectivity *const inBndConn, 
		                   const Field *const inState,
		        	   const Field *const inNode,
				   BoundaryConnectivity* outBndConn, 
				   Field* outState,
				   Field* outNode)
{
  LogToScreen(VERBOSE, "ShockFittingObj::processField()\n");  
}

//--------------------------------------------------------------------------//
  
void ShockFittingObj::createMeshData()
{
  MeshData::getInstance().createData <unsigned> ("NVT", 1);
  MeshData::getInstance().createData <unsigned> ("NBFACSH", 1);
  MeshData::getInstance().createData <unsigned> ("NFPOIN", 1);
  MeshData::getInstance().createData <unsigned> ("nPhanPoints",1);
  MeshData::getInstance().createData <unsigned> ("nBoundPhanPoints",1);
  MeshData::getInstance().createData <vector<string> > ("BoundariesMap",1);
  MeshData::getInstance().createData <vector<unsigned> > ("NPOIN", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NEDGE", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NELEM", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NBFAC", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NBPOIN", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NHOLE", 1);

  MeshData::getInstance().createData <vector<double> > ("CADDholes",1);

  MeshData::getInstance().createData <vector <int> >("NODCOD", 1);
  MeshData::getInstance().createData <vector <double> >("ZROE", 1);
  MeshData::getInstance().createData <vector <double> >("COOR", 1);
  MeshData::getInstance().createData <vector <int> >("BNDFAC", 1);
  MeshData::getInstance().createData <vector <int> >("CELNOD", 1);
  MeshData::getInstance().createData <vector <int> >("CELCEL", 1);
  MeshData::getInstance().createData <vector <int> >("EDGPTR", 1);
  MeshData::getInstance().createData <vector <int> >("NODPTR", 1);
  MeshData::getInstance().createData <vector <int> >("M12M0", 1);
  MeshData::getInstance().createData <vector <unsigned> >("M02M1", 1);
  MeshData::getInstance().createData <vector <int> >("ICLR", 1);

  MeshData::getInstance().createData <vector<double> >("MedianDualCellArea",1);
  MeshData::getInstance().createData <vector<unsigned> >("MedianDualCellNodes",1);
  MeshData::getInstance().createData <vector<unsigned> >("MedianDualCellNode",1);
  MeshData::getInstance().createData <vector<unsigned> >("MedianDualCellPtr",1);

  // store the primitive variables on the grid-points of the background
  // mesh for the current step and for the previous step
  MeshData::getInstance().createData <Array2D <double> >("primVariablesBkg",1);
  MeshData::getInstance().createData <Array2D <double> >("primVariablesBkgOld",1);

  MeshData::getInstance().createData <vector <int> >("NODCODBackup",1);
  MeshData::getInstance().createData <Array2D <int> >("NODPTRBackup",1);
  MeshData::getInstance().createData <Array2D <int> >("BNDFACBackup",1);

  MeshData::getInstance().createData <unsigned>("FirstRead",1);
  MeshData::getInstance().createData <stringstream>("FNAME");
  MeshData::getInstance().createData <string>("FNAMEBACK");

  MeshData::getInstance().createData <vector <double> >("firstResidual", 1);

  // vector used to store some CF connectivity data, the other
  // are the ones assigned to MeshData
  MeshData::getInstance().createData <vector <int> >("CF_boundaryNodes",1);
  MeshData::getInstance().createData <vector <string> >("CF_boundaryNames",1);
  MeshData::getInstance().createData <vector <int> >("CF_boundaryInfo",1);
  MeshData::getInstance().createData <vector <int> >("CF_boundaryPtr",1);
  MeshData::getInstance().createData <vector <int> >("CF_elementPtr",1);

  // store the connectivity, the mesh points state and coordinates
  // to exchange data with coofluid without using I/O
  MeshData::getInstance().createData <Connectivity> ("inConn");
  MeshData::getInstance().createData <BoundaryConnectivity> ("inbConn");
  MeshData::getInstance().createData <Field> ("inNodeField");
  MeshData::getInstance().createData <Field> ("inStateField");

  MeshData::getInstance().createData <Connectivity> ("outConn");
  MeshData::getInstance().createData <BoundaryConnectivity> ("outbConn");
  MeshData::getInstance().createData <Field> ("outNodeField");
  MeshData::getInstance().createData <Field> ("outStateField");

  MeshData::getInstance().setup();
}

//--------------------------------------------------------------------------//

void ShockFittingObj::createPhysicsData()
{
  PhysicsData::getInstance().createData <unsigned> ("NDOF", 1);

  PhysicsData::getInstance().createData <unsigned> ("NMOL",1);
  PhysicsData::getInstance().createData <unsigned> ("IE",1);
  PhysicsData::getInstance().createData <unsigned> ("IX",1);
  PhysicsData::getInstance().createData <unsigned> ("IY",1);
  PhysicsData::getInstance().createData <unsigned> ("IEV",1);

  PhysicsData::getInstance().createData <unsigned> ("NSP",1);
  PhysicsData::getInstance().createData 
                          <vector<string> > ("NAMESP",1);
  PhysicsData::getInstance().createData <vector<double> > ("MM",1);
  PhysicsData::getInstance().createData <vector<double> > ("HF",1);
  PhysicsData::getInstance().createData <vector<double> > ("THEV",1);
  PhysicsData::getInstance().createData <vector<double> > ("GAMS",1);
  PhysicsData::getInstance().createData 
                          <vector<string> > ("TYPEMOL",1);
  PhysicsData::getInstance().createData <vector<double> > ("RS",1);

  // these two values are set in ReferenceInfo object
  PhysicsData::getInstance().createData <double> ("GREF",1);
  PhysicsData::getInstance().createData <double> ("GM1REF",1);

  PhysicsData::getInstance().createData < unsigned > ("nShocks", 1);
  PhysicsData::getInstance().createData < unsigned > ("nSpecPoints", 1);
  PhysicsData::getInstance().createData 
                          <vector <unsigned> > ("nShockPoints", 1);
  PhysicsData::getInstance().createData 
			  <vector <unsigned> > ("nShockEdges", 1);
  PhysicsData::getInstance().createData 
			  <vector <string> > ("TypeSpecPoints", 1);
  PhysicsData::getInstance().createData <vector <string> > ("TYPESH", 1);
  PhysicsData::getInstance().createData <Array3D <double> > ("XYSH", 1);
  PhysicsData::getInstance().createData 
			  <Array3D <unsigned> > ("SHinSPPs", 1);
  PhysicsData::getInstance().createData 
			  <Array3D <double> > ("ZROESHuOLD", 1);
  PhysicsData::getInstance().createData 
			  <Array3D <double> > ("ZROESHdOLD", 1);
  PhysicsData::getInstance().createData <Array3D <double> > ("VSHNOR",1);
  PhysicsData::getInstance().createData <Array3D <double> > ("WSH",1);

  PhysicsData::getInstance().setup();
}

//--------------------------------------------------------------------------//

void ShockFittingObj::deleteMeshData()
{
  MeshData::getInstance().deleteData <unsigned> ("NVT");
  MeshData::getInstance().deleteData <unsigned> ("NBFACSH");
  MeshData::getInstance().deleteData <unsigned> ("NFPOIN");
  MeshData::getInstance().deleteData <unsigned> ("nPhanPoints");
  MeshData::getInstance().deleteData <unsigned> ("nBoundPhanPoints");
  MeshData::getInstance().deleteData <vector<string> > ("BoundariesMap");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NPOIN");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NEDGE");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NELEM");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NBFAC");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NBPOIN");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NHOLE");
  MeshData::getInstance().deleteData <vector<double> > ("CADDholes");

  MeshData::getInstance().deleteData <vector <int> >("NODCOD");
  MeshData::getInstance().deleteData <vector <double> >("ZROE");
  MeshData::getInstance().deleteData <vector <double> >("COOR");
  MeshData::getInstance().deleteData <vector <int> >("BNDFAC");
  MeshData::getInstance().deleteData <vector <int> >("CELNOD");
  MeshData::getInstance().deleteData <vector <int> >("CELCEL");
  MeshData::getInstance().deleteData <vector <int> >("EDGPTR");
  MeshData::getInstance().deleteData <vector <int> >("NODPTR");
  MeshData::getInstance().deleteData <vector <int> >("M12M0");
  MeshData::getInstance().deleteData <vector <unsigned> >("M02M1");
  MeshData::getInstance().deleteData <vector <int> >("ICLR");

  MeshData::getInstance().deleteData <vector<double> >("MedianDualCellArea");
  MeshData::getInstance().deleteData <vector<unsigned> >("MedianDualCellNodes");
  MeshData::getInstance().deleteData <vector<unsigned> >("MedianDualCellNode");
  MeshData::getInstance().deleteData <vector<unsigned> >("MedianDualCellPtr");

  MeshData::getInstance().deleteData <Array2D <double> >("primVariablesBkg"); 
  MeshData::getInstance().deleteData <Array2D <double> >("primVariablesBkgOld");

  MeshData::getInstance().deleteData <vector <int> >("NODCODBackup");
  MeshData::getInstance().deleteData <Array2D <int> >("NODPTRBackup");
  MeshData::getInstance().deleteData <Array2D <int> >("BNDFACBackup");

  MeshData::getInstance().deleteData <unsigned>("FirstRead");
  MeshData::getInstance().deleteData <stringstream >("FNAME");
  MeshData::getInstance().deleteData <string >("FNAMEBACK");

  MeshData::getInstance().deleteData <vector <double> >("firstResidual");

  MeshData::getInstance().deleteData <vector <int> >("CF_boundaryNodes");
  MeshData::getInstance().deleteData <vector <string> >("CF_boundaryNames");
  MeshData::getInstance().deleteData <vector <int> >("CF_boundaryInfo");
  MeshData::getInstance().deleteData <vector <int> >("CF_boundaryPtr");
  MeshData::getInstance().deleteData <vector <int> >("CF_elementPtr");

  MeshData::getInstance().deleteData <Connectivity> ("inConn");
  MeshData::getInstance().deleteData <BoundaryConnectivity> ("inbConn");
  MeshData::getInstance().deleteData <Field> ("inNodeField");
  MeshData::getInstance().deleteData <Field> ("inStateField");

  MeshData::getInstance().deleteData <Connectivity> ("outConn");
  MeshData::getInstance().deleteData <BoundaryConnectivity> ("outbConn");
  MeshData::getInstance().deleteData <Field> ("outNodeField");
  MeshData::getInstance().deleteData <Field> ("outStateField");

  MeshData::getInstance().unsetup();
}

//--------------------------------------------------------------------------//

void ShockFittingObj::deletePhysicsData()
{
  PhysicsData::getInstance().deleteData <unsigned> ("NDOF");

  PhysicsData::getInstance().deleteData <unsigned> ("NMOL");
  PhysicsData::getInstance().deleteData <unsigned> ("IE");
  PhysicsData::getInstance().deleteData <unsigned> ("IX");
  PhysicsData::getInstance().deleteData <unsigned> ("IY");
  PhysicsData::getInstance().deleteData <unsigned> ("IEV");

  PhysicsData::getInstance().deleteData <unsigned> ("NSP");
  PhysicsData::getInstance().deleteData <vector<string> > ("NAMESP");
  PhysicsData::getInstance().deleteData <vector<double> > ("MM");
  PhysicsData::getInstance().deleteData <vector<double> > ("HF");
  PhysicsData::getInstance().deleteData <vector<double> > ("THEV");
  PhysicsData::getInstance().deleteData <vector<double> > ("GAMS");
  PhysicsData::getInstance().deleteData <vector<string> > ("TYPEMOL");
  PhysicsData::getInstance().deleteData <vector<double> > ("RS");

  // these values are set in ReferenceInfo object
  PhysicsData::getInstance().deleteData <double> ("GREF");
  PhysicsData::getInstance().deleteData <double> ("GM1REF");

  PhysicsData::getInstance().deleteData < unsigned > ("nShocks");
  PhysicsData::getInstance().deleteData < unsigned > ("nSpecPoints");
  PhysicsData::getInstance().deleteData <vector <unsigned> > ("nShockPoints");
  PhysicsData::getInstance().deleteData <vector <unsigned> > ("nShockEdges");
  PhysicsData::getInstance().deleteData <vector <string> > ("TypeSpecPoints");
  PhysicsData::getInstance().deleteData <vector <string> > ("TYPESH");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("XYSH");
  PhysicsData::getInstance().deleteData <Array3D <unsigned> > ("SHinSPPs");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("ZROESHuOLD");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("ZROESHdOLD");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("VSHNOR");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("WSH");

  PhysicsData::getInstance().unsetup();
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
