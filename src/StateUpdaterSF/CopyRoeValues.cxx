// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/CopyRoeValues.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CopyRoeValues, StateUpdater> copyRoeValuesProv("CopyRoeValues");

//--------------------------------------------------------------------------//

CopyRoeValues::CopyRoeValues(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

CopyRoeValues::~CopyRoeValues()
{
  delete zroe1; delete zroe0;
}

//--------------------------------------------------------------------------//

void CopyRoeValues::setup()
{
  LogToScreen(VERBOSE, "CopyRoeValues::setup() => start\n");

  LogToScreen(VERBOSE, "CopyRoeValues::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues::unsetup()
{
  LogToScreen(VERBOSE, "CopyRoeValues::unsetup()\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues::update()
{
  LogToScreen(INFO, "CopyRoeValues::update()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned IA=0; IA<(*ndof); IA++) {
    (*zroe0)(IA,M12M0->at(IPOIN)) = (*zroe1)(IA,IPOIN);
   }
  }
}

//--------------------------------------------------------------------------//

void CopyRoeValues::setAddress()
{
  zroe0 = new Array2D<double>((*ndofmax),
                              (npoin->at(0)+2 * (*nshmax) * (*npshmax)),
                              &zroeVect->at(0));
  start = (*ndofmax) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  zroe1 = new Array2D<double>((*ndofmax),
                              (npoin->at(1)+2 * (*nshmax) * (*npshmax)),
                              &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void CopyRoeValues::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
}

//--------------------------------------------------------------------------//

void CopyRoeValues::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF"); 
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
