// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDpsPg.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MoveDpsPg, MoveDps> moveDpsPgProv("MoveDpsPg");

//----------------------------------------------------------------------------//

MoveDpsPg::MoveDpsPg(const std::string& objectName) :
 MoveDps(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDpsPg::~MoveDpsPg()
{
}

//----------------------------------------------------------------------------//

void MoveDpsPg::setup()
{
  LogToScreen(VERBOSE,"MoveDpsPg::setup() => start\n");

  LogToScreen(VERBOSE,"MoveDpsPg::setup() => end\n");
}

//----------------------------------------------------------------------------//

void MoveDpsPg::unsetup()
{
  LogToScreen(VERBOSE,"MoveDpsPg::unsetup()\n");
}

//----------------------------------------------------------------------------//

void MoveDpsPg::update()
{
  LogToScreen(INFO,"MoveDpsPg::update()\n");

  setMeshData();

  setPhysicsData();

  setAddress();

  logfile.Open(getClassName().c_str());

  // compute the dt max
  dt = 1e39;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock n. ", iShock, "\n");
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    unsigned iShockPoint = IV+1;
 
    ro = (*ZroeSh)(0,IV,ISH)*(*ZroeSh)(0,IV,ISH);
    help = pow((*ZroeSh)(2,IV,ISH),2) + pow((*ZroeSh)(3,IV,ISH),2);
    p = ((*gam)-1)/(*gam) * ((*ZroeSh)(0,IV,ISH) * (*ZroeSh)(1,IV,ISH)-0.5*help);
    a = sqrt((*gam)*p/ro);

    dum = 0;
    for(unsigned K=0; K<2; K++) { dum = dum + pow((*WSh)(K,IV,ISH),2); }
    WShMod = sqrt(dum);
    dum = 0;
    for(unsigned K=0; K<2; K++) { dum = dum + pow((*WSh)(K,IV,ISH),2); }
    WShMod = sqrt(dum);

    logfile("Shock point n. ", iShockPoint," speed: ", WShMod, "\n");

    dum = (*shrelax) * (*dxcell) * (*sndmin) /(a+WShMod);

    if(dt>dum) { dt = dum; }
   }
  }

  logfile("DT max: ", dt, "\n");

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock/Discontinuity n. ", iShock, "\n");
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    for(unsigned I=0; I<2; I++) {
     (*XYSh)(I,IV,ISH) = (*XYSh)(I,IV,ISH) + (*WSh)(I,IV,ISH) * dt;
    }
    logfile((*XYSh)(0,IV,ISH), " ", (*XYSh)(1,IV,ISH), "\n");
   }
  }
  logfile.Close();
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
