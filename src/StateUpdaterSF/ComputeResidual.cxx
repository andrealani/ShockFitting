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
  m_gasModel = "DummyModel"; // up-to-date Pg or TCneq 
  addOption("gasModel", &m_gasModel,
            "Gas model required for the conversion in primitive variables");
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

  m_paramToprimDimensional.ptr()->setup();

  LogToScreen(VERBOSE,"ComputeResidual::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ComputeResidual::unsetup()
{
  LogToScreen(VERBOSE,"ComputeResidual::unsetup()\n");

  m_normErr.ptr()->unsetup();

  m_paramToprimDimensional.ptr()->unsetup();
}

//--------------------------------------------------------------------------//

void ComputeResidual::configure(OptionMap& cmap, const std::string& prefix)
{
  StateUpdater::configure(cmap, prefix);

  // assign strings read on input.case file to the object computing the 
  // norm
  m_normErr.name() = "DiscrErrorNorm"+m_whichNorm;
  if(m_isItWeighted) { m_normErr.name() = m_normErr.name() + "weighted"; }

  // assign strings read in input.case to the object making the 
  // variable transformation
   m_paramToprimDimensional.name() = "Param2Prim"+m_gasModel+"Dimensional";

  if (ConfigFileReader::isFirstConfig()) {
   m_normErr.ptr().reset
     (SConfig::Factory<StateUpdater>::getInstance().
      getProvider(m_normErr.name())->create(m_normErr.name()));

   m_paramToprimDimensional.ptr().reset
     (SConfig::Factory<VariableTransformer>::getInstance().
      getProvider(m_paramToprimDimensional.name())->
      create(m_paramToprimDimensional.name()));
  }

  // configure the object computing the norm value
  configureDeps(cmap, m_normErr.ptr().get());

  // configure the object transforming the variables
  configureDeps(cmap, m_paramToprimDimensional.ptr().get());
}
 
//--------------------------------------------------------------------------//    

void ComputeResidual::update()
{
  LogToScreen(INFO,"ComputeResidual::update()\n");

  setMeshData();
  setPhysicsData();


  if(MeshData::getInstance().getIstep()==1) {

   // resize vector and arrays
   resizeArray();

   ofstream printNorm("SFconvergence.plt",ios::app);
   printNorm << "TITLE = Shock Fitting Convergence, norm: "<< m_whichNorm;
   if(m_isItWeighted) { printNorm << " weighted"; }
   printNorm << "\n";
   printNorm << "VARIABLES = ";
   printNorm << "Iter"  << " ";
   for(unsigned i=0; i<(*ndof); i++) { printNorm << "Res[" << i << "] "; }
   printNorm << endl;
   printNorm.close();
  }

  setAddress();

  // write the number of the current on the the output file
  if(MeshData::getInstance().getIstep()>1) {

   ofstream printNorm("SFconvergence.plt",ios::app);
   printNorm << MeshData::getInstance().getIstep() << " ";
   printNorm.close();
  }

  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   // assign the zroe of the IPOIN to the working vector m_zroe
   for(unsigned IV=0; IV<(*ndof); IV++) {
    m_zroe.at(IV) = (*zroe)(IV,IPOIN);
   }

   // convert the zroe variables in primitive variables
   m_paramToprimDimensional.ptr()->transform(&m_zroe,&m_XY,&m_prim);

   // assign the transformed prim variablesof the IPOIN-point  to 
   // the corrisponding array
   for(unsigned IV=0; IV<(*ndof); IV++) {        
    (*primBackgroundMesh)(IV,IPOIN) = m_prim.at(IV);
   }
  }

  if(MeshData::getInstance().getIstep()>1) {
   // compute the norm of the discretization error using the
   // primitive variables of the background mesh grid-points
   m_normErr.ptr()->update();
  }

  // update primBackgroundMeshOld with the values of the current step
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   for(unsigned K=0; K<(*ndof); K++) { 
    (*primBackgroundMeshOld)(K,IPOIN) = (*primBackgroundMesh)(K,IPOIN); 
   }
  }

  // de-allocate the dynamic arrays 
  freeArray();
}

//--------------------------------------------------------------------------//

void ComputeResidual::resizeArray()
{
  m_zroe.resize((*ndof));
  m_prim.resize((*ndof));
  m_XY.resize((*ndof),0); // it is not used
  primBackgroundMesh->resize((*ndof),npoin->at(0));
  primBackgroundMeshOld->resize((*ndof),npoin->at(0));
}

//--------------------------------------------------------------------------//

void ComputeResidual::setAddress()
{
  zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),
                               npoin->at(0),
                               &zroeVect->at(0));
}

//--------------------------------------------------------------------------//

void ComputeResidual::freeArray()
{
  delete zroe; 
}

//--------------------------------------------------------------------------//

void ComputeResidual::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  primBackgroundMesh = 
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkg");
  primBackgroundMeshOld = 
   MeshData::getInstance().getData <Array2D<double> >("primVariablesBkgOld");
}

//--------------------------------------------------------------------------//

void ComputeResidual::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace SHockFitting
