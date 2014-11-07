// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/CopyRoeValues1.hh"
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
ObjectProvider<CopyRoeValues1, CopyMaker> 
  copyRoeValuesProv("CopyRoeValues1");

//--------------------------------------------------------------------------//

CopyRoeValues1::CopyRoeValues1(const std::string& objectName) :
  CopyMaker(objectName)
{
}

//--------------------------------------------------------------------------//

CopyRoeValues1::~CopyRoeValues1()
{
  delete zroe1; delete zroe0;
}

//--------------------------------------------------------------------------//

void CopyRoeValues1::setup()
{
  LogToScreen(VERBOSE, "CopyRoeValues1::setup() => start\n");

  LogToScreen(VERBOSE, "CopyRoeValues1::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1::unsetup()
{
  LogToScreen(VERBOSE, "CopyRoeValues1::unsetup()\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1::copy()
{
  LogToScreen(INFO, "CopyRoeValues1::copy()\n");

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

void CopyRoeValues1::setAddress()
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

void CopyRoeValues1::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF"); 
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
