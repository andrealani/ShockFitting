// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <vector>
#include "MeshGeneratorSF/MeshInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Common/MeshData.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"


//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MeshInfo, MeshGenerator> meshInfoProv("MeshInfo");

//--------------------------------------------------------------------------//

MeshInfo::MeshInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{
  m_eps = 0;
  addOption("EPS",&m_eps,
	    "Distance between two shock faces");

  m_sndmin = 0;
  addOption("SNDMIN",&m_sndmin,
            "Max non dimensional distance of phantom nodes");

  m_dxcell = 0;
  addOption("DXCELL",&m_dxcell,
            "Length of the shock edges");

  m_ibak = 1;
  addOption("IBAK",&m_ibak,
            "Number of iterations before saving solution",true);

  m_naddholes = 0;
  addOption("Naddholes",&m_naddholes,
            "Number of hole points");

  m_caddholes = vector <double>();
  addOption("CADDholes",&m_caddholes,
            "Holes points coordinates");

  m_nproc = 1;
  addOption("NPROC",&m_nproc,
            "Number of processor");
}

//--------------------------------------------------------------------------//

MeshInfo::~MeshInfo()
{
}

//--------------------------------------------------------------------------//

void MeshInfo::setup()
{
  LogToScreen(VERBOSE, "MeshInfo::setup() => start\n");

  setMeshData();

  LogToScreen(VERBOSE, "MeshInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void MeshInfo::unsetup()
{
  LogToScreen(VERBOSE, "MeshInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void MeshInfo::generate()
{
  LogToScreen(INFO, "MeshInfo::generate()\n");

  *eps = m_eps;
  *sndmin = m_sndmin;
  *dxcell = m_dxcell;
  *ibak = m_ibak;
  *naddholes = m_naddholes;
  if ( (*naddholes) != 0 ) {
   caddholes->resize( 2 * (*naddholes) );
   for (unsigned i=0; i<caddholes->size(); i++) {
    caddholes->at(i) = m_caddholes.at(i);}
  }
  else {caddholes->resize(1); caddholes->at(0)=m_caddholes.at(0);}
  *nproc = m_nproc;
}

//--------------------------------------------------------------------------//

void MeshInfo::setMeshData()
{
  eps = MeshData::getInstance().getData <double> ("EPS");
  sndmin = MeshData::getInstance().getData <double> ("SNDMIN");
  dxcell = MeshData::getInstance().getData <double> ("DXCELL");
  ibak = MeshData::getInstance().getData <unsigned> ("IBAK");
  naddholes = MeshData::getInstance().getData <unsigned> ("Naddholes");
  caddholes = MeshData::getInstance().getData < std::vector<double> > ("CADDholes");
  nproc = MeshData::getInstance().getData <unsigned> ("NPROC");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting


