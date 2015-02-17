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
  setAddress();

  normValue.resize((*ndof),0);

  // L1w = (sum(i=1,N) |u^(n+1)-u^(n)|)/ (sum(i=1,N) |u^(1)-u^(0)|)
  // @param N number of mesh points
  
if(MeshData::getInstance().getIstep()>1) {
FILE* outfile = fopen("Zroe.txt","w");
for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    fprintf(outfile,"%32.16F %s",(*zroe)(K,IPOIN)," "); }
    fprintf(outfile,"%s","\n ");
   }
fclose(outfile);
outfile = fopen("ZroeOld.txt","w");
for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    fprintf(outfile,"%32.16F %s",(*zroeOld)(K,IPOIN)," "); }
    fprintf(outfile,"%s","\n ");
   }
fclose(outfile);
}


  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) {
    normValue.at(K) = normValue.at(K) + 
                      abs((*zroe)(K,IPOIN)-(*zroeOld)(K,IPOIN));
   }
  }

  // store the value of the first residual
  if(MeshData::getInstance().getIstep()==2) {
   firstResidualValue->resize((*ndof));
   for(unsigned K=0; K<(*ndof); K++) {
    firstResidualValue->at(K) = normValue.at(K) / npoin->at(1);
   }
  }

  // weight the residual value on the first residual
  for(unsigned K=0; K<(*ndof); K++) { 
    normValue.at(K) = normValue.at(K) / npoin->at(1) / firstResidualValue->at(K);
  }

  // define the fstream value printing the norm
  ofstream printNorm("SFconvergence.plt",ios::app);

  for(unsigned K=0; K<(*ndof); K++) {
   printNorm << normValue.at(K) << " "; }

  printNorm << endl;
  printNorm.close();

  // de-allocate the dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::setAddress()
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

void DiscrErrorNormL1weighted::freeArray()
{
  delete zroe; delete zroeOld;
}

//--------------------------------------------------------------------------//

void DiscrErrorNormL1weighted::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  zroeOldVect = MeshData::getInstance().getData <vector<double> >("ZROEOld");
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
