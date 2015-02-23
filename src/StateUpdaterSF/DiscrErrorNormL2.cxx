// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/DiscrErrorNormL2.hh"
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
ObjectProvider<DiscrErrorNormL2, StateUpdater> 
discrErrorNormL2Prov("DiscrErrorNormL2");

//--------------------------------------------------------------------------//

DiscrErrorNormL2::DiscrErrorNormL2(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

DiscrErrorNormL2::~DiscrErrorNormL2()
{
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::setup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL2::setup() => start\n");

  LogToScreen(VERBOSE,"DiscrErrorNormL2::setup() => end\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::unsetup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL2::unsetup()\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::update()
{
  LogToScreen(DEBUG_MIN,"DiscrErrorNormL2::update()\n");

  setMeshData();
  setPhysicsData();
  setAddress();

  unsigned currentNbMeshPoints;

  normValue.resize((*ndof),0);

  // L2 = ((sum(i=1,N) |u^(n+1)-u^(n)|^2)/ N)^1/2
  // @param N number of mesh points

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) + 
                      pow(abs((*zroe)(K,IPOIN)-(*zroeOld)(K,IPOIN)),2);
   }
  }

  // if the new number of points of the shocked mesh is smaller than
  // the old number of points of the shocked mesh, then the
  // residual will be computed on npoin->at(1)
  if(npoin->at(1)<(*npoinShockedMeshBkp)) {
   currentNbMeshPoints = npoin->at(1);
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     normValue.at(K) = normValue.at(K) +
                      abs((*zroe)(K,IPOIN)-(*zroeOld)(K,IPOIN));
    }
   }
  }

  // if the new number of points of the shocked mesh is larger than
  // the old number of points of the shocked mesh, then the
  // residual will be computed on npoinShockedMeshBkp
  else {
   currentNbMeshPoints = (*npoinShockedMeshBkp);
   for(unsigned IPOIN=0; IPOIN<(*npoinShockedMeshBkp); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     normValue.at(K) = normValue.at(K) +
                      abs((*zroe)(K,IPOIN)-(*zroeOld)(K,IPOIN));
    }
   }
  }
        
  // compute the norm
  for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = sqrt(normValue.at(K)/currentNbMeshPoints); 
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

void DiscrErrorNormL2::setAddress()
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

void DiscrErrorNormL2::freeArray()
{
  delete zroe; delete zroeOld;
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::setMeshData()
{
  npoinShockedMeshBkp =
    MeshData::getInstance().getData <unsigned>("NPOINshockedMeshBkp");   
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  zroeOldVect = MeshData::getInstance().getData <vector<double> >("ZROEOld");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
