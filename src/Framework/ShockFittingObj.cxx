// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ShockFittingObj.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
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
            "List of the names of mesh generator files");

  m_cNormalVector.name() = "CoNorm";
  addOption("CoNorm", &m_cNormalVector,
            "Object computing normal vector");

  m_fRemeshing = vector<PAIR_TYPE(Remeshing)>();
  addOption("RemeshingList",&m_fRemeshing,
	    "List of the names of field remeshing");

  m_wMesh = vector<PAIR_TYPE(WritingMesh)>();
  addOption("WritingMeshList",&m_wMesh,
            "List of the names of objects writing");

  m_fConverter = vector<PAIR_TYPE(Converter)>();
  addOption("ConverterList",&m_fConverter,
            "List of the names of converter objects");

  m_CFDSolver.name() = "COOLFluiD";

  m_sUpdater = vector<PAIR_TYPE(StateUpdater)>();
  addOption("StateUpdaterList",&m_sUpdater,
            "List of the names of state updaters");

  m_cState.name() = "ComputeStateDps";
  addOption("ComputeStateDps", &m_cState,
            "Object updating solution by enforcing RH relations");

  m_moveDps.name() = "MoveDps";
  addOption("MoveDps", &m_moveDps,
            "Object moving shock points");
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

    // create the remeshing
    createList<Remeshing>(m_fRemeshing);

    // create the normal vector computing
    m_cNormalVector.ptr().reset(SConfig::Factory<CoNorm>::getInstance().
                          getProvider(m_cNormalVector.name())
                          ->create(m_cNormalVector.name()));

    // create the writing objects
    createList<WritingMesh>(m_wMesh);

    // create the converter objects
    createList<Converter>(m_fConverter);

    // create the object calling COOLFluiD
    m_CFDSolver.ptr().reset(SConfig::Factory<CFDSolver>::getInstance().
                          getProvider(m_CFDSolver.name())
                          ->create(m_CFDSolver.name()));

    // create the state updaters list
    createList<StateUpdater>(m_sUpdater);

    // create the solution updating
    m_cState.ptr().reset(SConfig::Factory<ComputeStateDps>::getInstance().
                          getProvider(m_cState.name())
                          ->create(m_cState.name()));

    // create shock moving shock points
    m_moveDps.ptr().reset(SConfig::Factory<MoveDps>::getInstance().
                          getProvider(m_moveDps.name())
                          ->create(m_moveDps.name()));
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

  // configure the file mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    configureDeps (cmap, m_mGenerator[i].ptr().get());
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    configureDeps (cmap, m_fRemeshing[i].ptr().get());
  }

  // configure the normal vector computing
  configureDeps (cmap, m_cNormalVector.ptr().get());

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
   configureDeps (cmap, m_wMesh[i].ptr().get());
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
   configureDeps (cmap, m_fConverter[i].ptr().get());
  }

  // configure the object calling COOLFluiD
  configureDeps (cmap,m_CFDSolver.ptr().get());

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
   configureDeps (cmap, m_sUpdater[i].ptr().get());
  }

  // configure the object updating solution
  configureDeps (cmap,m_cState.ptr().get());

  // configure the object moving shock points
  configureDeps (cmap,m_moveDps.ptr().get());

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
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->setup();
  } 

  // configure the file mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    m_mGenerator[i].ptr()->setup();
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    m_fRemeshing[i].ptr()->setup();
  }

  // configure the normal vector computing
  m_cNormalVector.ptr()->setup();

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
    m_wMesh[i].ptr()->setup();
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
    m_fConverter[i].ptr()->setup();
  }

  // configure the object calling COOLFluiD
  m_CFDSolver.ptr()->setup();

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
    m_sUpdater[i].ptr()->setup();
  }

  // configure the computing state
  m_cState.ptr()->setup();

  // configure object moving shock points
  m_moveDps.ptr()->setup();

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

  // configure the file mesh generator
  for (unsigned i = 0; i < m_mGenerator.size(); ++i) {
    m_mGenerator[i].ptr()->unsetup();
  }

  // configure the remeshing
  for (unsigned i = 0; i < m_fRemeshing.size(); ++i) {
    m_fRemeshing[i].ptr()->unsetup();
  }

  // configure the normal vector computing
  m_cNormalVector.ptr()->unsetup();

  // configure the writing objects
  for (unsigned i = 0; i < m_wMesh.size(); ++i) {
    m_wMesh[i].ptr()->unsetup();
  }

  // configure the converter objects
  for (unsigned i = 0; i < m_fConverter.size(); ++i) {
    m_fConverter[i].ptr()->unsetup();
  }

  // configure the object calling COOLFluiD
  m_CFDSolver.ptr()->unsetup();

  // configure the state updating
  for (unsigned i = 0; i < m_sUpdater.size(); ++i) {
    m_sUpdater[i].ptr()->unsetup();
  }

  // configure the object updating solution
  m_cState.ptr()->unsetup();

  // configure object moving shock points
  m_moveDps.ptr()->setup();

  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::process()
{
  LogToScreen(VERBOSE, "ShockFittingObj::process()\n");
}

//--------------------------------------------------------------------------//

void ShockFittingObj::createMeshData()
{
  MeshData::getInstance().createData <unsigned> ("NVT", 1);
  MeshData::getInstance().createData <unsigned> ("NBFACSH", 1);
  MeshData::getInstance().createData <unsigned> ("NFPOIN", 1);
  MeshData::getInstance().createData <unsigned> ("nPhanPoints",1);
  MeshData::getInstance().createData <unsigned> ("nBoundPhanPoints",1);
  MeshData::getInstance().createData <vector<unsigned> > ("NPOIN", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NEDGE", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NELEM", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NBFAC", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NBPOIN", 2);
  MeshData::getInstance().createData <vector<unsigned> > ("NHOLE", 1);

  MeshData::getInstance().createData <double> ("EPS",1);
  MeshData::getInstance().createData <double> ("SNDMIN",1);
  MeshData::getInstance().createData <double> ("DXCELL",1);
  MeshData::getInstance().createData <double> ("SHRELAX",1);
  MeshData::getInstance().createData <unsigned> ("IBAK",1);
  MeshData::getInstance().createData <unsigned> ("Naddholes",1);
  MeshData::getInstance().createData <vector<double> > ("CADDholes",1);
  MeshData::getInstance().createData <unsigned> ("NPROC",1);

  MeshData::getInstance().createData <vector <int> >("NODCOD", 1);
  MeshData::getInstance().createData <vector <double> >("ZROE", 1);
  MeshData::getInstance().createData <vector <double> >("COOR", 1);
  MeshData::getInstance().createData <vector <int> >("BNDFAC", 1);
  MeshData::getInstance().createData <vector <int> >("CELNOD", 1);
  MeshData::getInstance().createData <vector <int> >("CELCEL", 1);
  MeshData::getInstance().createData <vector <int> >("EDGPTR", 1);
  MeshData::getInstance().createData <vector <int> >("NODPTR", 1);
  MeshData::getInstance().createData <vector <int> >("M12M0", 1);
  MeshData::getInstance().createData <vector <int> >("M02M1", 1);

  MeshData::getInstance().createData <unsigned>("FirstRead",1);
  MeshData::getInstance().createData <vector<string> >("FNAME",1);
  MeshData::getInstance().createData <vector<string> >("FNAMEBACK",1);

  MeshData::getInstance().setup();
}

//--------------------------------------------------------------------------//

void ShockFittingObj::createPhysicsData()
{
  
  PhysicsData::getInstance().createData <unsigned> ("NDIM", 1);
  // these two values are read by PhysicsInfo object
  PhysicsData::getInstance().createData <double> ("GAM", 1);
  PhysicsData::getInstance().createData <double> ("GM1", 1);
  PhysicsData::getInstance().createData <unsigned> ("NDOF", 1);
  PhysicsData::getInstance().createData <unsigned> ("NDOFMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NSHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NPSHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NSPMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NESHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NADDHOLESMAX", 1);
  PhysicsData::getInstance().createData 
                          <vector<string> > ("MODEL", 1);
  PhysicsData::getInstance().createData 
                          <vector<string> > ("MIXTURE", 1);
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

  // these two values are read by ReferenceInfo object
  PhysicsData::getInstance().createData <double> ("RgasFreeStream",1);
  PhysicsData::getInstance().createData <double> ("GamFreeStream",1);
  PhysicsData::getInstance().createData <double> ("PREF",1);
  PhysicsData::getInstance().createData <double> ("TREF",1);
  PhysicsData::getInstance().createData <double> ("UREF",1);
  PhysicsData::getInstance().createData <double> ("RHOREF",1);
  PhysicsData::getInstance().createData <double> ("LREF",1);
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
  MeshData::getInstance().deleteData <vector<unsigned> > ("NPOIN");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NEDGE");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NELEM");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NBFAC");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NBPOIN");
  MeshData::getInstance().deleteData <vector<unsigned> > ("NHOLE");

  MeshData::getInstance().deleteData <double> ("EPS");
  MeshData::getInstance().deleteData <double> ("SNDMIN");
  MeshData::getInstance().deleteData <double> ("DXCELL");
  MeshData::getInstance().deleteData <double> ("SHRELAX");
  MeshData::getInstance().deleteData <unsigned> ("IBAK");
  MeshData::getInstance().deleteData <unsigned> ("Naddholes");
  MeshData::getInstance().deleteData <vector<double> > ("CADDholes");
  MeshData::getInstance().deleteData <unsigned> ("NPROC");

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

  MeshData::getInstance().deleteData <unsigned>("FirstRead");
  MeshData::getInstance().deleteData <vector<string> >("FNAME");
  MeshData::getInstance().deleteData <vector<string> >("FNAMEBACK");

  MeshData::getInstance().unsetup();
}

//--------------------------------------------------------------------------//

void ShockFittingObj::deletePhysicsData()
{
  PhysicsData::getInstance().unsetup();
  
  PhysicsData::getInstance().deleteData <unsigned> ("NDIM");
  PhysicsData::getInstance().deleteData <unsigned> ("NDOF");
  PhysicsData::getInstance().deleteData <unsigned> ("NDOFMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NSHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NPSHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NSPMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NESHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NADDHOLESMAX");
  // these values are read by PhysicsInfo object
  PhysicsData::getInstance().deleteData <double> ("GAM");
  PhysicsData::getInstance().deleteData <double> ("GM1");

  PhysicsData::getInstance().deleteData <vector<string> > ("MODEL");
  PhysicsData::getInstance().deleteData <vector<string> > ("MIXTURE");

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

  // these values are read by ReferenceInfo object
  PhysicsData::getInstance().deleteData <double> ("RgasFreeStream");
  PhysicsData::getInstance().deleteData <double> ("GamFreeStream");
  PhysicsData::getInstance().deleteData <double> ("PREF");
  PhysicsData::getInstance().deleteData <double> ("TREF");
  PhysicsData::getInstance().deleteData <double> ("UREF");
  PhysicsData::getInstance().deleteData <double> ("RHOREF");
  PhysicsData::getInstance().deleteData <double> ("LREF");
  PhysicsData::getInstance().deleteData <double> ("GREF");
  PhysicsData::getInstance().deleteData <double> ("GM1REF");

  PhysicsData::getInstance().deleteData < unsigned > ("nShocks");
  PhysicsData::getInstance().deleteData < unsigned > ("nSpecPoints");
  PhysicsData::getInstance().deleteData <vector <unsigned> > ("nShockPoints");
  PhysicsData::getInstance().deleteData <vector <unsigned> > ("nShockEdges");
  PhysicsData::getInstance().deleteData <vector <string> > ("TypeSpecPoints");
  PhysicsData::getInstance().deleteData <vector <string> > ("TYPESH");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("ZROESHuOLD");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("ZROESHdOLD");
  PhysicsData::getInstance().deleteData <Array3D <unsigned> > ("SHinSPPs");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("VSHNOR");
  PhysicsData::getInstance().deleteData <Array3D <double> > ("WSH");

  PhysicsData::getInstance().unsetup();
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
