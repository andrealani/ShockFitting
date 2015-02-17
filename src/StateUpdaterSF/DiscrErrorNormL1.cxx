// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium

// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/DiscrErrorNormL1.hh"
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
ObjectProvider<DiscrErrorNormL1, StateUpdater> 
discrErrorNormL1Prov("DiscrErrorNormL1");

//--------------------------------------------------------------------------//

DiscrErrorNormL1::DiscrErrorNormL1(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

DiscrErrorNormL1::~DiscrErrorNormL1()
{
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL1::setup() => start\n");

  LogToScreen(VERBOSE,"DiscrErrorNormL1::setup() => end\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::unsetup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL1::unsetup()\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::update()
{
  LogToScreen(DEBUG_MIN,"DiscrErrorNormL1::update()\n");

  setMeshData();
  setPhysicsData();
  setAddress();

  normValue.resize((*ndof),0);

  // L1 = (sum(i=1,N) |u^(n+1)-u^(n)|)/ N
  // @param N number of mesh points

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) + 
                      abs((*zroe)(K,IPOIN)-(*zroeOld)(K,IPOIN));
   }
  }

  for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) / npoin->at(1); 
  }

  // define the fstream value printing the norm
  ofstream printNorm("SFconvergence.plt",ios::app);

  for(unsigned K=0; K<(*ndof); K++) {
   printNorm << normValue.at(K) << " "; }

  printNorm << endl;
  printNorm.close();

  // de-allocate the dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setAddress()
{
  unsigned start = PhysicsInfo::getnbDofMax() *
                   (npoin->at(0) + 2 *
                    PhysicsInfo::getnbShMax() *
                    PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),
                              (npoin->at(1) + 2 *
                               PhysicsInfo::getnbShMax() *
                               PhysicsInfo::getnbShPointsMax()),
                               &zroeVect->at(start));
  zroeOld = new Array2D <double>  (PhysicsInfo::getnbDofMax(),
                                  (npoin->at(1) + 2 *
                                   PhysicsInfo::getnbShMax() *
                                   PhysicsInfo::getnbShPointsMax()),
                                   &zroeOldVect->at(0));
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::freeArray()
{
  delete zroe; delete zroeOld;
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  zroeOldVect = MeshData::getInstance().getData <vector<double> >("ZROEOld");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
