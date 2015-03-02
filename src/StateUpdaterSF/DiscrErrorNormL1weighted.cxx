// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium

// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/DiscrErrorNormL1weighted.hh"
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
ObjectProvider<DiscrErrorNormL1weighted, StateUpdater> 
discrErrorNormL1weightedProv("DiscrErrorNormL1weighted");

//--------------------------------------------------------------------------//

DiscrErrorNormL1weighted::DiscrErrorNormL1weighted(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

DiscrErrorNormL1weighted::~DiscrErrorNormL1weighted()
{
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::setup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL1weighted::setup() => start\n");

  LogToScreen(VERBOSE,"DiscrErrorNormL1weighted::setup() => end\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::unsetup()
{
  LogToScreen(VERBOSE,"DiscrErrorNormL1weighted::unsetup()\n");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::update()
{
  LogToScreen(DEBUG_MIN,"DiscrErrorNormL1weighted::update()\n");

  setMeshData();
  setPhysicsData();

  normValue.resize((*ndof),0);

  // L1w = (sum(i=1,N) |u_(n+1)-u_(n)|)/ (sum(i=1,N) |u_(1)-u_(0)|)
  // @param N number of mesh points
  
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) + 
                      abs((*primBackgroundMesh)(K,IPOIN)-
                          (*primBackgroundMeshOld)(K,IPOIN));
   }
  }

  // store the value of the first residual
  if(MeshData::getInstance().getIstep()==2) {
   firstResidualValue->resize((*ndof));
   for(unsigned K=0; K<(*ndof); K++) {
    firstResidualValue->at(K) = normValue.at(K) / npoin->at(0);
   }
  }

  // weight the residual value on the first residual
  for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) / npoin->at(0) / firstResidualValue->at(K);
  }

  // define the fstream value printing the norm
  ofstream printNorm("SFconvergence.plt",ios::app);

  for(unsigned K=0; K<(*ndof); K++) {
   printNorm << normValue.at(K) << " "; }

  printNorm << endl;
  printNorm.close();
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  primBackgroundMesh =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkg");
  primBackgroundMeshOld =
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkgOld");
  firstResidualValue = 
   MeshData::getInstance().getData <vector<double> >("firstResidual");
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
