// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ShockFittingObj.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/ConfigFileReader.hh"

#include "Common/MeshData.hh"
#include "Common/PhysicsData.hh"
#include "MathTools/Array2D.hh"

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
  m_vTransformer = vector<PAIR_TYPE(VariableTransformer)>();
  addOption("VariableTransformerList",&m_vTransformer, 
	    "List of the names of variable transformers"); 
  
  m_fInterpolator = vector<PAIR_TYPE(FieldInterpolator)>();
  addOption("FieldInterpolatorList",&m_fInterpolator, 
	    "List of the names of field interpolators"); 
  
  m_fProcessing = vector<PAIR_TYPE(FileProcessing)>();
  addOption("FileProcessingList",&m_fProcessing, 
	    "List of the names of file processing"); 

  m_fMeshGenerator = vector<PAIR_TYPE(MeshGenerator)>();
  addOption("MeshGeneratorList",&m_fMeshGenerator,
            "List of the names of mesh generator files");
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
    
    // create the variable transformers
    createList<VariableTransformer>(m_vTransformer);
    
    // create the field interpolators
    createList<FieldInterpolator>(m_fInterpolator);
    
    // create the file processing
    createList<FileProcessing>(m_fProcessing);

    // create the file mesh generator
    createList<MeshGenerator>(m_fMeshGenerator);
  }
  
  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    configureDeps (cmap, m_vTransformer[i].ptr().get());
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    configureDeps (cmap, m_fInterpolator[i].ptr().get());
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    configureDeps (cmap, m_fProcessing[i].ptr().get());
  }
 
  // configure the file mesh generator
  for (unsigned i = 0; i < m_fMeshGenerator.size(); ++i) {
    configureDeps (cmap, m_fMeshGenerator[i].ptr().get());
  }
  
  LogToScreen(VERBOSE, "ShockFittingObj::configure() => end\n");
}

//--------------------------------------------------------------------------//

void ShockFittingObj::setup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::setup() => start\n");
 
  createMeshData();

  createPhysicsData();

  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    m_vTransformer[i].ptr()->setup();
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->setup();
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->setup();
  } 

  // configure the file mesh generator
  for (unsigned i = 0; i < m_fMeshGenerator.size(); ++i) {
    m_fMeshGenerator[i].ptr()->setup();
  }
  
  LogToScreen(VERBOSE, "ShockFittingObj::setup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::unsetup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => start\n");

  deleteMeshData();

  deletePhysicsData();
  
  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    m_vTransformer[i].ptr()->unsetup();
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->unsetup();
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->unsetup();
  }

  // configure the file mesh generator
  for (unsigned i = 0; i < m_fMeshGenerator.size(); ++i) {
    m_fMeshGenerator[i].ptr()->unsetup();
  }
  
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
  MeshData::getInstance().createData <unsigned> ("NPOIN", 1);
  MeshData::getInstance().createData <unsigned> ("NBFAC", 1);
  MeshData::getInstance().createData <unsigned> ("NBPOIN", 1);
  MeshData::getInstance().createData <unsigned> ("NFPOIN", 1);
  MeshData::getInstance().createData <unsigned> ("NHOLE", 1);

  MeshData::getInstance().createData <double> ("EPS",1);
  MeshData::getInstance().createData <double> ("SNDMIN",1);
  MeshData::getInstance().createData <double> ("DXCELL",1);
  MeshData::getInstance().createData <unsigned> ("IBAK",1);
  MeshData::getInstance().createData <unsigned> ("Naddholes",1);
  MeshData::getInstance().createData < std::vector<double> > ("CADDholes",1);
  MeshData::getInstance().createData <unsigned> ("NPROC",1);

  MeshData::getInstance().createData <vector <int> >("NODCOD", 1);
  MeshData::getInstance().createData <vector <double> >("ZROE", 1);
  MeshData::getInstance().createData <vector <double> >("COOR", 1);
  MeshData::getInstance().createData <Array2D <int> >("BNDFAC", 1);
  MeshData::getInstance().createData <Array2D <int> >("CELNOD", 1);
  MeshData::getInstance().createData <Array2D <int> >("CELCEL", 1);
  MeshData::getInstance().createData <Array2D <int> >("EDGPTR", 1);
  MeshData::getInstance().createData <Array2D <int> >("NODPTR", 1);
}

//--------------------------------------------------------------------------//

void ShockFittingObj::createPhysicsData()
{
  PhysicsData::getInstance().createData <unsigned> ("NDIM", 1);
  PhysicsData::getInstance().createData <unsigned> ("NDOF", 1);
  PhysicsData::getInstance().createData <unsigned> ("NDOFMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NSHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NPSHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NSPMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NESHMAX", 1);
  PhysicsData::getInstance().createData <unsigned> ("NADDHOLESMAX", 1);
  PhysicsData::getInstance().createData <char> ("MODEL", 1);
  PhysicsData::getInstance().createData <char> ("MIXTURE", 1);

}

//--------------------------------------------------------------------------//

void ShockFittingObj::deleteMeshData()
{
  MeshData::getInstance().deleteData <unsigned> ("NPOIN");
  MeshData::getInstance().deleteData <unsigned> ("NBFAC");
  MeshData::getInstance().deleteData <unsigned> ("NBPOIN");
  MeshData::getInstance().deleteData <unsigned> ("NFPOIN");
  MeshData::getInstance().deleteData <unsigned> ("NHOLE");

  MeshData::getInstance().deleteData <double> ("EPS");
  MeshData::getInstance().deleteData <double> ("SNDMIN");
  MeshData::getInstance().deleteData <double> ("DXCELL");
  MeshData::getInstance().deleteData <unsigned> ("IBAK");
  MeshData::getInstance().deleteData <unsigned> ("Naddholes");
  MeshData::getInstance().deleteData < std::vector<double> > ("CADDholes");
  MeshData::getInstance().deleteData <unsigned> ("NPROC");

  MeshData::getInstance().deleteData <vector <int> >("NODCOD");
  MeshData::getInstance().deleteData <vector <double> >("ZROE");
  MeshData::getInstance().deleteData <vector <double> >("COOR");
  MeshData::getInstance().deleteData <Array2D <int> >("BNDFAC");
  MeshData::getInstance().deleteData <Array2D <int> >("CELNOD");
  MeshData::getInstance().deleteData <Array2D <int> >("CELCEL");
  MeshData::getInstance().deleteData <Array2D <int> >("EDGPTR");
  MeshData::getInstance().deleteData <Array2D <int> >("NODPTR");
}

//--------------------------------------------------------------------------//

void ShockFittingObj::deletePhysicsData()
{
  PhysicsData::getInstance().deleteData <unsigned> ("NDIM");
  PhysicsData::getInstance().deleteData <unsigned> ("NDOF");
  PhysicsData::getInstance().deleteData <unsigned> ("NDOFMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NSHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NPSHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NSPMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NESHMAX");
  PhysicsData::getInstance().deleteData <unsigned> ("NADDHOLESMAX");
  PhysicsData::getInstance().deleteData <char> ("MODEL");
  PhysicsData::getInstance().deleteData <char> ("MIXTURE");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
