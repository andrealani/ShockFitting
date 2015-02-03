// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/MeshRestoring.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MeshRestoring, CopyMaker>
  meshRestoringProv("MeshRestoring");

//--------------------------------------------------------------------------//

MeshRestoring::MeshRestoring(const std::string& objectName) :
  CopyMaker(objectName)
{
}

//--------------------------------------------------------------------------//

MeshRestoring::~MeshRestoring()
{
}

//--------------------------------------------------------------------------//

void MeshRestoring::setup()
{
  LogToScreen(VERBOSE,"MeshRestoring::setup() => start\n");

  LogToScreen(VERBOSE,"MeshRestoring::setup() => end\n");
}

//--------------------------------------------------------------------------//

void MeshRestoring::unsetup()
{
  LogToScreen(VERBOSE,"MeshRestoring::unsetup()\n");
}

//--------------------------------------------------------------------------//

void MeshRestoring::copy()
{
  LogToScreen(INFO,"MeshRestoring::copy()\n");

  setMeshData();

  setAddress();

  // restore nbpoin
  nbpoin->at(0) = nbpoin->at(2);

  // restore nbfac
  nbfac->at(0) = nbfac->at(2);

  // restore bndfac
  for(unsigned IFACE=0; IFACE<nbfac->at(2); IFACE++) {
   for(unsigned I=0; I<3; I++) {
    (*bndfac)(I,IFACE) = (*bndfacBackup)(I,IFACE);
   }
  }

  // restore nodptr
  for(unsigned IBPOIN=0; IBPOIN<nbpoin->at(2); IBPOIN++) {
   for(unsigned I=0; I<3; I++) {
    (*nodptr)(IBPOIN,I) = (*nodptrBackup)(IBPOIN,I);
   }
  }

  // restore nodcod
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   nodcod->at(IPOIN) = nodcodBackup->at(IPOIN);
  }

  // de-allocate the dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void MeshRestoring::setAddress()
{
  unsigned totsize = nbfac->at(0) + 2 *
                     PhysicsInfo::getnbShMax() *
                     PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(0));
  nodptr = new Array2D<int> (nbpoin->at(0),3, &nodptrVect->at(0));
}

//--------------------------------------------------------------------------//

void MeshRestoring::freeArray()
{
  delete bndfac; delete nodptr;
}

//--------------------------------------------------------------------------//

void MeshRestoring::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbpoin = MeshData::getInstance().getData <vector<unsigned> > ("NBPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  nodptrVect = MeshData::getInstance().getData <vector<int> >("NODPTR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  nodcodBackup =
     MeshData::getInstance().getData <vector<int> >("NODCODBackup");
  nodptrBackup =
     MeshData::getInstance().getData <Array2D<int> >("NODPTRBackup");
  bndfacBackup =
     MeshData::getInstance().getData <Array2D<int> >("BNDFACBackup");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
