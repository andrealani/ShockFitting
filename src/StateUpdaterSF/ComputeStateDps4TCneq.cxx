// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/ComputeStateDps4TCneq.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
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
ObjectProvider<ComputeStateDps4TCneq, ComputeStateDps>
 computeStateDpsTCneqProv("ComputeStateDps4TCneq");

//----------------------------------------------------------------------------//

ComputeStateDps4TCneq::ComputeStateDps4TCneq(const std::string& objectName) :
  ComputeStateDps(objectName)
{
  TotnbShockPoints = 0;
}

//----------------------------------------------------------------------------//

ComputeStateDps4TCneq::~ComputeStateDps4TCneq()
{
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::setup()
{
  LogToScreen(VERBOSE,"ComputeStateDps4TCneq::setup() => start\n");

  LogToScreen(VERBOSE,"ComputeStateDps4TCneq::setup() => end\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::unsetup()
{
  LogToScreen(VERBOSE,"ComputeStateDps4TCneq::unsetup()\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::update()
{
  LogToScreen(INFO,"ComputeStateDps4TCneq::update()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setDiscSpeedSize();

  unsigned I;

  // create object of CoShock class
  CoShock computenewStateForShock;

  // create object of CoDc class
  CoDc computenewStateForDc;

  alphad.resize( (*nsp) );
  alphau.resize( (*nsp) );
  xd.resize(4);
  xu.resize(4);

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {

   unsigned ivalue = ISH+1;
   logfile("Shock/Disc. n. ",ivalue, "\n");

   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    ++TotnbShockPoints;

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
     computenewStateForShock.callCoShock(xd,xu,R2);
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

    evd = evu;

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

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::recoverDownState(unsigned IV, unsigned ISH)
{
  double ZRHOSH = 0.0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   ZRHOSH = ZRHOSH + (*ZroeShd)(ISP,IV,ISH); // sqrt(rho)
  }
  // density
  xd.at(0) = ZRHOSH*ZRHOSH;

  double RHOHF = 0.0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   RHOHF = RHOHF + (*ZroeShd)(ISP,IV,ISH)* hf->at(ISP);
  }
  // formation enthalpy
  RHOHF = RHOHF*ZRHOSH;  // (rho*Hf)

  evd = (*ZroeShd)((*IEV),IV,ISH)/ZRHOSH;

  xd.at(3) = ( -(*ZroeShd)((*IX),IV,ISH) * dy + (*ZroeShd)((*IY),IV,ISH) * dx );
  xd.at(3) = xd.at(3)/ZRHOSH; // tangential
  xd.at(2) = ( (*ZroeShd)((*IX),IV,ISH) * dx + (*ZroeShd)((*IY),IV,ISH) * dy );
  xd.at(2) = xd.at(2)/ZRHOSH; // normal

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   // chemical species concentrations
   alphad.at(ISP) = (*ZroeShd)(ISP,IV,ISH)/ZRHOSH;
  }

  help = pow( (*ZroeShd)((*IX),IV,ISH) , 2) + pow( (*ZroeShd)((*IY),IV,ISH) , 2);
  // pressure
  xd.at(1) = ((*gref)-1)/(*gref) * (ZRHOSH *  (*ZroeShd)((*IE),IV,ISH) - 0.50*help
             - RHOHF - ZRHOSH* (*ZroeShd)((*IEV),IV,ISH));

  R2 = sqrt((*gref)*xd.at(1)/xd.at(0)) + 0.50 * ((*gref)-1.0) * xd.at(2);
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::recoverUpState(unsigned IV, unsigned ISH)
{ 
  double ZRHOSH = 0.0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   ZRHOSH = ZRHOSH + (*ZroeShu)(ISP,IV,ISH); // sqrt(rho)
  }
  // density
  xu.at(0) = ZRHOSH*ZRHOSH;
  
  double RHOHF = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) { 
   RHOHF = RHOHF + (*ZroeShu)(ISP,IV,ISH)* hf->at(ISP);
  }
  // formation enthalpy 
  RHOHF = RHOHF*ZRHOSH;  // (rho*Hf)
  
  evu = (*ZroeShu)((*IEV),IV,ISH)/ZRHOSH;
  
  xu.at(3) = ( -(*ZroeShu)((*IX),IV,ISH) * dy + (*ZroeShu)((*IY),IV,ISH) * dx );
  xu.at(3) = xu.at(3)/ZRHOSH; // tangential
  xu.at(2) = ( (*ZroeShu)((*IX),IV,ISH) * dx + (*ZroeShu)((*IY),IV,ISH) * dy );
  xu.at(2) = xu.at(2)/ZRHOSH; // normal
  
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   // chemical species concentrations
   alphau.at(ISP) = (*ZroeShu)(ISP,IV,ISH)/ZRHOSH;
  }
  
  help = pow( (*ZroeShu)((*IX),IV,ISH) , 2) + pow( (*ZroeShu)((*IY),IV,ISH) , 2);
  // pressure
  xu.at(1) = ((*gref)-1.0)/(*gref) * (ZRHOSH *  (*ZroeShu)((*IE),IV,ISH) - 0.50*help
             - RHOHF - ZRHOSH* (*ZroeShu)((*IEV),IV,ISH));
} 

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::computeDownState(unsigned IV, unsigned ISH)
{
  vector <double> zsv(*nsp);
  double z1, z2, z3, z4, z5, kine, enthalpy;

  z1 = sqrt(xd.at(0));
  kine = 0.5 * (xd.at(2)*xd.at(2)+xd.at(3)*xd.at(3));
  double hfv = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   hfv = hfv + alphau.at(ISP) * hf->at(ISP);
  }
  enthalpy = (*gref)/((*gref)-1) * xd.at(1)/xd.at(0) + kine + hfv + evd;
  z2 = z1*enthalpy;
  z3 = z1 * (xd.at(2) * dx - xd.at(3) * dy);
  z4 = z1 * (xd.at(2) * dy + xd.at(3) * dx);
  z5 = z1 * evd;

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   zsv.at(ISP) = alphau.at(ISP)*z1;
  }

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   (*ZroeShd)(ISP,IV,ISH) = zsv.at(ISP);
  }

  (*ZroeShd)((*IE),IV,ISH) = z2;
  (*ZroeShd)((*IX),IV,ISH) = z3;
  (*ZroeShd)((*IY),IV,ISH) = z4;
  (*ZroeShd)((*IEV),IV,ISH) = z5;
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::computeUpState(unsigned IV, unsigned ISH)
{ 
  vector <double> zsv(*nsp);
  double z1, z2, z3, z4, z5, kine, enthalpy;
  
  z1 = sqrt(xu.at(0));
  kine = 0.5 * (pow(xu.at(2),2)+pow(xu.at(3),2));
  double hfv = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   hfv = hfv + alphau.at(ISP) * hf->at(ISP);
  }
  enthalpy = (*gref)/((*gref)-1) * xu.at(1)/xu.at(0) + kine + hfv + evd;
  z2 = z1*enthalpy;
  z3 = z1 * (xu.at(2) * dx - xu.at(3) * dy);
  z4 = z1 * (xu.at(2) * dy + xu.at(3) * dx);
  z5 = z1 * evu;
  
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   zsv.at(ISP) = alphau.at(ISP)*z1;
  }
  
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   (*ZroeShu)(ISP,IV,ISH) = zsv.at(ISP);
  }

  (*ZroeShu)((*IE),IV,ISH) = z2;
  (*ZroeShu)((*IX),IV,ISH) = z3;
  (*ZroeShu)((*IY),IV,ISH) = z4;
  (*ZroeShu)((*IEV),IV,ISH) = z5;
}

//----------------------------------------------------------------------------//

void ComputeStateDps4TCneq::saveDownState(unsigned IV, unsigned ISH)
{
  for(unsigned I=0; I<(*ndof); I++) {
   (*ZroeShdOld)(I,IV,ISH) = (*ZroeShd)(I,IV,ISH);
  }
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

