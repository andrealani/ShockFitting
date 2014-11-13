// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/CopyRoeValues2.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CopyRoeValues2, CopyMaker> 
 copyRoeValues2Prov("CopyRoeValues2");

//--------------------------------------------------------------------------//

CopyRoeValues2::CopyRoeValues2(const std::string& objectName) :
  CopyMaker(objectName)
{
}

//--------------------------------------------------------------------------//

CopyRoeValues2::~CopyRoeValues2()
{
}

//--------------------------------------------------------------------------//

void CopyRoeValues2::setup()
{
  LogToScreen(VERBOSE, "CopyRoeValues2::setup() => start\n");

  LogToScreen(VERBOSE, "CopyRoeValues2::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues2::unsetup()
{
  LogToScreen(VERBOSE, "CopyRoeValues2::unsetup()\n");
}

//--------------------------------------------------------------------------//

void CopyRoeValues2::copy()
{
  LogToScreen(INFO, "CopyRoeValues2::copy()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  unsigned IB;
  unsigned ILIST = 2 * PhysicsInfo::getnbShMax() * 
                       PhysicsInfo::getnbShPointsMax();
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

void CopyRoeValues2::setAddress()
{
  zroe0 = new Array2D<double>((*ndof),
                              (npoin->at(0) + 2 * 
                               PhysicsInfo::getnbShPointsMax() * 
                               PhysicsInfo::getnbShMax()),
                              &zroeVect->at(0));
  start = (*ndof) * (npoin->at(0) + 2 * 
           PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax());
  zroe1 = new Array2D<double>((*ndof),
                              (npoin->at(1)+ 2 *
                               PhysicsInfo::getnbShPointsMax() *
                               PhysicsInfo::getnbShMax()),
                              &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void CopyRoeValues2::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
}

//--------------------------------------------------------------------------//

void CopyRoeValues2::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

