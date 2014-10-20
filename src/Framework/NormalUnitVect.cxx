// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NormalUnitVect.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "RemeshingSF/CoNorm.hh"
#include "RemeshingSF/CoNorm4B.hh"
#include "RemeshingSF/CoNormPG.hh"
#include "RemeshingSF/CoNorm4Ar.hh"
#include "RemeshingSF/CoNorm4TCneq.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<NormalUnitVect, Remeshing> normalUnitVectProv("NormalUnitVect");

//--------------------------------------------------------------------------//

NormalUnitVect::NormalUnitVect(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

NormalUnitVect::~NormalUnitVect()
{
}

//--------------------------------------------------------------------------//

void NormalUnitVect::remesh()
{
  LogToScreen(INFO,"NormalUnitVect::remesh()\n");

  setPhysicsData();

  CoNorm* computeNormal = new CoNorm("CoNorm");

  if(getModel()=="B"&&(*ndof)==1) {
   CoNorm4B computeNormalVector4B("CoNorm4B");
   computeNormal=&computeNormalVector4B;

   computeNormal->setup();
   computeNormal->remesh();
   computeNormal->unsetup();
  }

   else if (getModel()=="PG" && (*ndof)==4) {
   CoNormPG computeNormalVectorPG("CoNormPG");
   computeNormal=&computeNormalVectorPG;

   computeNormal->setup();
   computeNormal->remesh();
   computeNormal->unsetup();
  }

  else if (getModel()=="Cneq" && getMixture()=="ar4") {
   CoNorm4Ar computeNormalVector4Ar("CoNorm4Ar");
   computeNormal=&computeNormalVector4Ar;

   computeNormal->setup();
   computeNormal->remesh();
   computeNormal->unsetup();
  }

  else if (getModel()=="TCneq") {
   CoNorm4TCneq computeNormalVector4TCneq("CoNorm4TCneq");
   computeNormal=&computeNormalVector4TCneq;

   computeNormal->setup();
   computeNormal->remesh();
   computeNormal->unsetup();
  }
}

//--------------------------------------------------------------------------//

string NormalUnitVect::getModel() const {return model->at(0);}

//--------------------------------------------------------------------------//

string NormalUnitVect::getMixture() const {return mixture->at(0);}

//--------------------------------------------------------------------------//

unsigned NormalUnitVect::getNbDof() const {return *ndof;}

//--------------------------------------------------------------------------//

void NormalUnitVect::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  model = PhysicsData::getInstance().getData <vector<string> > ("MODEL"); 
  mixture = PhysicsData::getInstance().getData <vector<string> > ("MIXTURE");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
