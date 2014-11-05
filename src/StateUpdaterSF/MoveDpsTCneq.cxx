// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDpsTCneq.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MoveDpsTCneq, MoveDps> moveDpsTCneqProv("MoveDpsTCneq");

//----------------------------------------------------------------------------//

MoveDpsTCneq::MoveDpsTCneq(const std::string& objectName) :
 MoveDps(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDpsTCneq::~MoveDpsTCneq()
{
}

//----------------------------------------------------------------------------//

void MoveDpsTCneq::setup()
{
  LogToScreen(VERBOSE,"MoveDpsTCneq::setup() => start\n");

  LogToScreen(VERBOSE,"MoveDpsTCneq::setup() => end\n");
}

//----------------------------------------------------------------------------//

void MoveDpsTCneq::unsetup()
{
  LogToScreen(VERBOSE,"MoveDpsTCneq::unsetup()\n");
}

//----------------------------------------------------------------------------//

void MoveDpsTCneq::update()
{
  LogToScreen(INFO,"MoveDpsTCneq::update()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  double ZRHO;
  double RHOHF;

  logfile.Open(getClassName().c_str());

  // compute the dt max
  dt = 1e39;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock n. ", iShock, "\n");
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    unsigned iShockPoint = IV+1;
    ZRHO = 0; RHOHF = 0;
    for(unsigned ISP=0; ISP<(*nsp); ISP++) {
     ZRHO = ZRHO + (*ZroeSh)(ISP,IV,ISH);
     RHOHF = RHOHF + (*ZroeSh)(ISP,IV,ISH) * hf->at(ISP);
    }
    RHOHF = RHOHF * ZRHO; 

    ro = ZRHO*ZRHO;
    help = pow((*ZroeSh)((*ix),IV,ISH),2)+pow((*ZroeSh)((*iy),IV,ISH),2);
    p = ((*gref)-1.0)/(*gref) * (ZRHO * (*ZroeSh)((*ie),IV,ISH) - 
         0.50 * help - RHOHF - (*ZroeSh)((*iev),IV,ISH) * ZRHO);
    a = sqrt((*gref)*p/ro);

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
