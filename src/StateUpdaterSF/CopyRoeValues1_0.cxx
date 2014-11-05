// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/CopyRoeValues1_0.hh"
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
ObjectProvider<CopyRoeValues1_0, StateUpdater> 
 copyRoeValues1_0Prov("CopyRoeValues1_0");

//--------------------------------------------------------------------------//

CopyRoeValues1_0::CopyRoeValues1_0(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

CopyRoeValues1_0::~CopyRoeValues1_0()
{
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::setup()
{
  LogToScreen(VERBOSE, "CopyRoeValues1_0::setup() => start\n");

  LogToScreen(VERBOSE, "CopyRoeValues1_0::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::unsetup()
{
  LogToScreen(VERBOSE, "CopyRoeValues1_0::unsetup()\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::update()
{
  LogToScreen(INFO, "CopyRoeValues1_0::update()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  unsigned IB;
  unsigned ILIST = 2 * (*nshmax) * (*npshmax);
  unsigned m_npoin0 = npoin->at(0)+1;

  for(unsigned IPOIN=0; IPOIN<ILIST; IPOIN++) {
   IB = M02M1->at(m_npoin0+IPOIN);
   if(IB != 0) {
    for(unsigned IA=0; IA<(*ndof); IA++) {
     (*zroe1)(IA,IB-1) = (*zroe0)(IA,IPOIN+npoin->at(0));
    }
   }
  }
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::setAddress()
{
  zroe0 = new Array2D<double>((*ndof),
                              (npoin->at(0) + 2 * (*nshmax) * (*npshmax)),
                              &zroeVect->at(0));
  start = (*ndof) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  zroe1 = new Array2D<double>((*ndof),
                              (npoin->at(1)+ 2 * (*nshmax) * (*npshmax)),
                              &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
}

//--------------------------------------------------------------------------//

void CopyRoeValues1_0::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

