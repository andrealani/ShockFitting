// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/ComputeResidual.hh"
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
ObjectProvider<ComputeResidual, StateUpdater> computeResidualProv("ComputeResidual");

//--------------------------------------------------------------------------//

ComputeResidual::ComputeResidual(const std::string& objectName) :
  StateUpdater(objectName)
{
  m_normErr.name() = "DummyNorm";
  m_whichNorm = "L1orL2";
  addOption("wichNorm", &m_whichNorm,
            "Specifies which ype of norm will be used");
  m_isItWeighted = false;
  addOption("isItWeighted", &m_isItWeighted,
            "Specifies if the used norm is weighted");
}

//--------------------------------------------------------------------------//

ComputeResidual::~ComputeResidual()
{
}

//--------------------------------------------------------------------------//

void ComputeResidual::setup()
{
  LogToScreen(VERBOSE,"ComputeResidual::setup() => start\n");

  m_normErr.ptr()->setup();

  LogToScreen(VERBOSE,"ComputeResidual::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ComputeResidual::unsetup()
{
  LogToScreen(VERBOSE,"ComputeResidual::unsetup()\n");

  m_normErr.ptr()->unsetup();
}

//--------------------------------------------------------------------------//

void ComputeResidual::configure(OptionMap& cmap, const std::string& prefix)
{
  StateUpdater::configure(cmap, prefix);

  // assign strings read on input.case file to the object computing the 
  // norm
  m_normErr.name() = "DiscrErrorNorm"+m_whichNorm;
  if(m_isItWeighted) { m_normErr.name() = m_normErr.name() + "weighted"; }

  if (ConfigFileReader::isFirstConfig()) {
   m_normErr.ptr().reset(SConfig::Factory<StateUpdater>::getInstance().
                             getProvider(m_normErr.name())
                             ->create(m_normErr.name()));
  }

  // configure the object computing the norm value
  configureDeps(cmap, m_normErr.ptr().get());
}
 
//--------------------------------------------------------------------------//    

void ComputeResidual::update()
{
  LogToScreen(INFO,"ComputeResidual::update()\n");

  setMeshData();
  setPhysicsData();

  if(MeshData::getInstance().getIstep()==1) {
   // resize zroeOldVect if the first step
   unsigned totsize = npoin->at(1) +
             2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();
   zroeOldVect->resize(PhysicsInfo::getnbDofMax() * totsize);

   ofstream printNorm("SFconvergence.plt",ios::app);
   printNorm << "TITLE = Convergence of Shock Fitting norm: "<< m_whichNorm;
   if(m_isItWeighted) { printNorm << " weighted"; }
   printNorm << "\n";
   printNorm << "VARIABLES = ";
   printNorm << "Iter"  << " ";
   for(unsigned i=0; i<(*ndof); i++) { printNorm << "Res[" << i << "] "; }
   printNorm << endl;
   printNorm.close();
  }

  setAddress();

  if(MeshData::getInstance().getIstep()>1) {
   ofstream printNorm("SFconvergence.plt",ios::app);
   printNorm << MeshData::getInstance().getIstep() << " ";
   printNorm.close();
  }

  if(MeshData::getInstance().getIstep()>1) {
   // compute the norm of the discretization error
   m_normErr.ptr()->update();
  }

  // update zroeOld with the values of the current step
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) { 
    (*zroeOld)(K,IPOIN) = (*zroe)(K,IPOIN); 
   }
  }

  // de-allocate the dynamic arrays 
  freeArray();
}

//--------------------------------------------------------------------------//

void ComputeResidual::setAddress()
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

void ComputeResidual::freeArray()
{
  delete zroe; delete zroeOld;
}

//--------------------------------------------------------------------------//

void ComputeResidual::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  zroeOldVect = MeshData::getInstance().getData <vector<double> >("ZROEOld");
}

//--------------------------------------------------------------------------//

void ComputeResidual::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace SHockFitting
