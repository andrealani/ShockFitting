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

  normValue.resize((*ndof),0);
  for(unsigned K=0; K<(*ndof); K++) { normValue.at(K) = 0; }

  // L2 = ((sum(i=1,N) |u_(n+1)-u_(n)|^2)/ N)^1/2
  // @param N number of mesh points

  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) + 
                      pow(abs((*primBackgroundMesh)(K,IPOIN)-
                              (*primBackgroundMeshOld)(K,IPOIN)),2);
   }
  }

  // compute the norm
  for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = sqrt(normValue.at(K)/npoin->at(0)); 
  }

  string fileConv = MeshData::getInstance().getResultsDir() + "/SFconvergence.plt";

  // define the fstream value printing the norm
  ofstream printNorm(fileConv.c_str(),ios::app);
  
  for(unsigned K=0; K<(*ndof); K++) {
   printNorm << normValue.at(K) << " "; }

  printNorm << endl;
  printNorm.close();
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  primBackgroundMesh =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkg");
  primBackgroundMeshOld =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkgOld");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL2::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
