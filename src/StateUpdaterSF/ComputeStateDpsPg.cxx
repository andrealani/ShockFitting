// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/ComputeStateDpsPg.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"
#include "StateUpdaterSF/CoDc.hh"
#include "StateUpdaterSF/CoShock.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ComputeStateDpsPg, ComputeStateDps>
 computeStateDpsPgProv("ComputeStateDpsPg");

//----------------------------------------------------------------------------//

ComputeStateDpsPg::ComputeStateDpsPg(const std::string& objectName) :
 ComputeStateDps(objectName)
{
}

//----------------------------------------------------------------------------//

ComputeStateDpsPg::~ComputeStateDpsPg()
{
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::setup()
{
  LogToScreen(VERBOSE,"ComputeStateDpsPg::setup() => start\n");

  LogToScreen(VERBOSE,"ComputeStateDpsPg::setup() => end\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::unsetup()
{
  LogToScreen(VERBOSE,"ComputeStateDpsPg::unsetup()\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::update()
{
  LogToScreen(INFO,"ComputeStateDpsPg::update()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setDiscSpeedSize();

  // create object of CoShock class
  CoShock computenewStateForShock;

  // create object of CoDc class
  CoDc computenewStateForDc;

  xd.resize(4);
  xu.resize(4);

  unsigned I;

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {

   unsigned ivalue = ISH+1;
   logfile("Shock/Disc. n. ",ivalue, "\n");

   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    ++TotnbShockPoints;

    R2.resize(nShockPoints->at(ISH),(*nShocks));

    I=IV;
    dx = (*vShNor)(0,I,ISH);
    dy = (*vShNor)(1,I,ISH);

    // upload downstream status
    recoverDownState(IV,ISH);

    // upload upstream status
    recoverUpState(IV,ISH);

    // initialize discontinuity speed
    WS = 0.0;

   if(typeSh->at(ISH)=="S") {
    computenewStateForShock.callCoShock(xd,xu,R2(IV,ISH));
    xd = computenewStateForShock.getnewDownValues();
    WS = computenewStateForShock.getnewDiscSpeed();
   }

   if(typeSh->at(ISH)=="D") {
    computenewStateForDc.callCoDc(xd,xu);
    xd = computenewStateForDc.getnewDownValues();
    xu = computenewStateForDc.getnewUpValues();
    WS = computenewStateForDc.getnewDiscSpeed();
   }

   // enforce tangential component equality for the shock case
   if(typeSh->at(ISH)=="S") { xd.at(3) = xu.at(3); }

   // save old downstream status
   saveDownState(IV, ISH);

   // compute downstream variables and assing the new values
   // to Zroe array (for the shock case)
   computeDownState(IV, ISH);

   // compute upstream variables and assing the new values 
   // to Zroe array (for the shock case)
   computeUpState(IV, ISH);

   // set the new discontinuity speed
   (*WSh)(0,IV,ISH) = WS * dx;
   (*WSh)(1,IV,ISH) = WS * dy;

   ivalue = IV+1;
   logfile("S/D point nr. ", ivalue, " Speed: ");
   logfile((*WSh)(0,IV,ISH), ", ", (*WSh)(1,IV,ISH), "\n");
   }
  }
  logfile.Close();
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::recoverDownState(unsigned IV, unsigned ISH)
{
  xd.at(3) = -(*ZroeShd)(2,IV,ISH) * dy + (*ZroeShd)(3,IV,ISH) * dx;
  xd.at(3) = xd.at(3)/(*ZroeShd)(0,IV,ISH); // tangential
  xd.at(2) = (*ZroeShd)(2,IV,ISH) * dx + (*ZroeShd)(3,IV,ISH) * dy;
  xd.at(2) = xd.at(2)/(*ZroeShd)(0,IV,ISH); // normal
  xd.at(0) = (*ZroeShd)(0,IV,ISH)*(*ZroeShd)(0,IV,ISH); // density
  help = pow((*ZroeShd)(2,IV,ISH),2) + pow((*ZroeShd)(3,IV,ISH),2);
  xd.at(1) = ((*gref)-1)/(*gref) * ((*ZroeShd)(0,IV,ISH)
              * (*ZroeShd)(1,IV,ISH)) - 0.5 * help; // pressure
  R2(IV,ISH) = sqrt((*gref)*xd.at(1)/xd.at(0)) + 0.5 * ((*gref)-1) * xd.at(2);

  logfile("Zd(1) ",(*ZroeShd)(0,IV,ISH),"\n");
  logfile("Zd(2) ",(*ZroeShd)(1,IV,ISH),"\n");
  logfile("Zd(3) ",(*ZroeShd)(2,IV,ISH),"\n");
  logfile("Zd(4) ",(*ZroeShd)(3,IV,ISH),"\n");
}


//----------------------------------------------------------------------------//

void ComputeStateDpsPg::recoverUpState(unsigned IV, unsigned ISH)
{
  xu.at(3) = -(*ZroeShu)(2,IV,ISH) * dy + (*ZroeShu)(3,IV,ISH) * dx;
  xu.at(3) = xu.at(3)/(*ZroeShu)(0,IV,ISH); // tangential
  xu.at(2) = (*ZroeShu)(2,IV,ISH) * dx + (*ZroeShu)(3,IV,ISH) * dy;
  xu.at(2) = xu.at(2)/(*ZroeShu)(0,IV,ISH); // normal
  xu.at(0) = (*ZroeShu)(0,IV,ISH)*(*ZroeShu)(0,IV,ISH); // density
  help = pow((*ZroeShu)(2,IV,ISH),2) + pow((*ZroeShu)(3,IV,ISH),2);
  xu.at(1) = ((*gref)-1)/(*gref) * ((*ZroeShu)(0,IV,ISH)
              * (*ZroeShu)(1,IV,ISH)) - 0.5 * help; // pressure

  logfile("Zu(1) ",(*ZroeShu)(0,IV,ISH),"\n");
  logfile("Zu(2) ",(*ZroeShu)(1,IV,ISH),"\n");
  logfile("Zu(3) ",(*ZroeShu)(2,IV,ISH),"\n");
  logfile("Zu(4) ",(*ZroeShu)(3,IV,ISH),"\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::computeDownState(unsigned IV, unsigned ISH)
{
  double z1, z2, z3, z4, hh, kine;

  z1 = sqrt(xd.at(0));
  kine = 0.5 * (pow(xd.at(2),2)+pow(xd.at(3),2));
  hh = (*gref)/((*gref)-1) * xd.at(1)/xd.at(0) + kine;
  z2 = z1*hh;
  z3 = z1 * (xd.at(2) * dx - xd.at(3) * dy);
  z4 = z1 * (xd.at(2) * dy + xd.at(3) * dx);

  (*ZroeShd)(0,IV,ISH) = z1;
  (*ZroeShd)(1,IV,ISH) = z2;
  (*ZroeShd)(2,IV,ISH) = z3;
  (*ZroeShd)(3,IV,ISH) = z4;
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::computeUpState(unsigned IV, unsigned ISH)
{
  double z1, z2, z3, z4, hh, kine;

  z1 = sqrt(xu.at(0));
  kine = 0.5 * (pow(xu.at(2),2)+pow(xu.at(3),2));
  hh = (*gref)/((*gref)-1) * xu.at(1)/xu.at(0) + kine;
  z2 = z1*hh;
  z3 = z1 * (xu.at(2) * dx - xu.at(3) * dy);
  z4 = z1 * (xu.at(2) * dy + xu.at(3) * dx);

  (*ZroeShu)(0,IV,ISH) = z1;
  (*ZroeShu)(1,IV,ISH) = z2;
  (*ZroeShu)(2,IV,ISH) = z3;
  (*ZroeShu)(3,IV,ISH) = z4;
}

//----------------------------------------------------------------------------//

void ComputeStateDpsPg::saveDownState(unsigned IV, unsigned ISH)
{
  for(unsigned I=0; I<(*ndof); I++) {
   (*ZroeShdOld)(I,IV,ISH) = (*ZroeShd)(I,IV,ISH);
  }
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
