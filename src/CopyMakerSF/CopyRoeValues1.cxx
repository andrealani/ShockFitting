// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/CopyRoeValues1.hh"
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

/*
ifstream var;
stringstream pathvar;
pathvar.str(string());
if(MeshData::getInstance().getIstep()<10){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step0000"<<MeshData::getInstance().getIstep()<<"/Var/dcopy1.var";
}
else if (MeshData::getInstance().getIstep()>=10 &&
         MeshData::getInstance().getIstep()<100){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step000"<<MeshData::getInstance().getIstep()<<"/Var/dcopy1.var";
}


string path = pathvar.str();
var.open(path.c_str());
if(var.fail()) { cout << "Step000" << MeshData::getInstance().getIstep() << "Failed opening dcopy1.var" << endl;
}


for(unsigned I=0;I<npoin->at(1);I++){
for(unsigned IA=0;IA<(*ndof);IA++){
var >> (*zroe1)(IA,I);}}
var.close();
*/

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned IA=0; IA<(*ndof); IA++) {
    // M12M0 has filled with indeces that start from 1
    // M12M0(1:2*NSHMAX*NPSHMAX)
    (*zroe0)(IA,M12M0->at(IPOIN+1)-1) = (*zroe1)(IA,IPOIN);
   }
  }
}

//--------------------------------------------------------------------------//

void CopyRoeValues1::setAddress()
{
  zroe0 = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                              (npoin->at(0)+ 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                              &zroeVect->at(0));
  start = PhysicsInfo::getnbDofMax() * 
          (npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax());
  zroe1 = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                              (npoin->at(1)+ 2 *
                              PhysicsInfo::getnbShMax() * 
                              PhysicsInfo::getnbShPointsMax()),
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
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
