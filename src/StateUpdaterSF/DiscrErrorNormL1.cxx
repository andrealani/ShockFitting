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

  normValue.resize((*ndof),0);
  for(unsigned K=0; K<(*ndof); K++) { normValue.at(K) = 0; }

  // L1 = (sum(i=1,N) |u_(n+1)-u_(n)|)/ N
  // @param N number of mesh points
  
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) +
                      abs((*primBackgroundMesh)(K,IPOIN)-
                          (*primBackgroundMeshOld)(K,IPOIN));
   }
  }

  // compute the residual
  for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) / npoin->at(0); 
  }

  // define the fstream value printing the norm
  ofstream printNorm("SFconvergence.plt",ios::app);

  for(unsigned K=0; K<(*ndof); K++) {
   printNorm << normValue.at(K) << " "; }

  printNorm << endl;
  printNorm.close();
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  primBackgroundMesh =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkg");
  primBackgroundMeshOld =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkgOld");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
