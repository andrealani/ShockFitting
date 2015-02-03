// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/FixStateSps.hh"
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
ObjectProvider<FixStateSps, StateUpdater> fixStateSpsProv("FixStateSps");

//--------------------------------------------------------------------------//

FixStateSps::FixStateSps(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

FixStateSps::~FixStateSps()
{
}

//--------------------------------------------------------------------------//

void FixStateSps::setup()
{
  LogToScreen(VERBOSE,"FixStateSps::setup() => start\n");

  LogToScreen(VERBOSE,"FixStateSps::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FixStateSps::unsetup()
{
  LogToScreen(VERBOSE,"FixStateSps::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FixStateSps::update()
{
  LogToScreen(INFO,"FixStateSps::update()\n");

  logfile.Open(getClassName().c_str());

  setMeshData();
  setPhysicsData();

  setAddress();

  vector<double> varZ((*ndof));
  vector<double> avarZ((*ndof));

  logfile("Z variables on entry in FixStateSps\n");
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock n. ",iShock, "\n");
   for(unsigned I=0; I<nShockPoints->at(ISH); I++) {
    unsigned iSpecPoint = I+1;
    logfile(iSpecPoint, " ");
    for(unsigned IV=0; IV<(*ndof); IV++) { logfile((*ZroeShu)(IV,I,ISH), " "); }
    for(unsigned IV=0; IV<(*ndof); IV++) { logfile((*ZroeShd)(IV,I,ISH), " "); }
    logfile("\n");
   }
  }

  for(unsigned ISPPNTS=0; ISPPNTS<(*nSpecPoints); ISPPNTS++) {
   
   if       (typeSpecPoints->at(ISPPNTS) == "IPX" ||
             typeSpecPoints->at(ISPPNTS) == "IPY" )
    // boundary special point: inlet point (floating along X or Y)
    { fixIPXandIPYspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "OPX" ||
             typeSpecPoints->at(ISPPNTS) == "OPY" )
    // boundary special point: outlet point (floating along X or Y)
    { fixOPXandOPYspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "WPNRX" ||
             typeSpecPoints->at(ISPPNTS) == "WPNRY" )
    // boundary special point: wall point without reflection 
    // (floating along X or Y)
    { fixWPNRXandWPNRYspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "TP" )
    // boundary special point: triple point
    { fixTPspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "QP" )
    // boundary special point: quadruple point
    { fixQPspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "RRX" )
    // boundary special point: quadruple point
    { fixRRXspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "EP" )
    // boundary special point: end point in a supersonic zone
    { fixEPspecPoints(ISPPNTS); }

   else if  (typeSpecPoints->at(ISPPNTS) == "C" )
    // boundary special point: connection between two shocks
    { fixCspecPoints(ISPPNTS); }

   else { cout << "FixStateSps::error => Condition not implemented\n";
          exit(1); }
  }

  // compute changes in the states of the shock points
  WWS = 0;
  for(unsigned IV=0; IV<(*ndof); IV++) {
   varZ .at(IV) = 0;
   avarZ.at(IV) = 0;
  }

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Z variations in the downstream zone of the shock\n");
   logfile("Shock n. ",iShock, "\n");
   for(unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for(unsigned IV=0; IV<(*ndof); IV++) {
     varZ.at(IV) = (*ZroeShd)(IV,I,ISH)-(*ZroeShdOld)(IV,I,ISH);
     avarZ.at(IV) = avarZ.at(IV) + sqrt(pow(varZ.at(IV),2));
    }
    WWS = WWS + sqrt(pow((*WSh)(0,I,ISH),2) + pow((*WSh)(1,I,ISH),2));

    logfile(I, " ");
    for(unsigned IV=0; IV<(*ndof); IV++) { logfile(varZ.at(IV), " "); }
    logfile("\n");
   }
  }

  logfile("Average Z values:\n");
  for(unsigned IV=0; IV<(*ndof); IV++) { logfile(avarZ.at(IV), " "); }
  logfile("\n");

  // save old upstream states and assign computes upstream state
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned IP=0; IP<nShockPoints->at(ISH); IP++) {
    for(unsigned IV=0; IV<(*ndof); IV++) {
     (*ZroeShdOld)(IV,IP,ISH) = (*ZroeShd)(IV,IP,ISH);
     (*ZroeShuOld)(IV,IP,ISH) = (*ZroeShu)(IV,IP,ISH);
    }
   }
  }

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void FixStateSps::fixIPXandIPYspecPoints(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  // restore downstream and upstream states
  for(unsigned I=0; I<(*ndof); I++) {
   (*ZroeShu)(I,IP.at(0),ISH.at(0)) = (*ZroeShuOld)(I,IP.at(0),ISH.at(0));
  }

  // set the shock speed to 0
  for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) 
   { (*WSh)(I,IP.at(0),ISH.at(0))=0; }
}

//--------------------------------------------------------------------------//

void FixStateSps::fixOPXandOPYspecPoints(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  // set the shock speed
  WS = sqrt(pow((*WSh)(0,IP.at(0),ISH.at(0)),2) + 
            pow((*WSh)(1,IP.at(0),ISH.at(0)),2));
  dx = (*WSh)(0,IP.at(0),ISH.at(0))/WS;
  dy = (*WSh)(1,IP.at(0),ISH.at(0))/WS;

  if      (typeSpecPoints->at(ISPPNTS) == "OPX") { 
   (*WSh)(0,IP.at(0),ISH.at(0)) = WS/dx;
   (*WSh)(1,IP.at(0),ISH.at(0)) = 0;             }
  else if (typeSpecPoints->at(ISPPNTS) == "OPY") {
   (*WSh)(0,IP.at(0),ISH.at(0)) = 0;
   (*WSh)(1,IP.at(0),ISH.at(0)) = WS/dy;       
  }  
}

//--------------------------------------------------------------------------//

void FixStateSps::fixWPNRXandWPNRYspecPoints(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  // set the shock speed
  WS = sqrt(pow((*WSh)(0,IP.at(0),ISH.at(0)),2) + 
            pow((*WSh)(1,IP.at(0),ISH.at(0)),2));
  dx = (*WSh)(0,IP.at(0),ISH.at(0))/WS;
  dy = (*WSh)(1,IP.at(0),ISH.at(0))/WS;
  
  if      (typeSpecPoints->at(ISPPNTS) == "WPNRX") { 
   (*WSh)(0,IP.at(0),ISH.at(0)) = WS/dx;
   (*WSh)(1,IP.at(0),ISH.at(0)) = 0;               }
  else if (typeSpecPoints->at(ISPPNTS) == "WPNRY") { 
   (*WSh)(0,IP.at(0),ISH.at(0)) = 0;
   (*WSh)(1,IP.at(0),ISH.at(0)) = WS/dy;           }   
}

//--------------------------------------------------------------------------//

void FixStateSps::fixTPspecPoints(unsigned ISPPNTS)
{
/*  // create CoUTP object
  CoUTP computeTPstates;

  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  // ISH.at(2) mach stem
  // ISH.at(3) contact discontinuity
  setShockIndeces(4,ISPPNTS);

  // set zone states

  // state (1): state between incident shock and mach stem
  setState("Up",1);

  // state (2): state between incident shock and reflected shock
  setState("Up",2);

  // state (3): state between reflected shock and contact discontinuity
  // set R23 Riemann invariant also
  setState("Down",3);

  // state (4): state between contact discontinuity and mach stem
  // set R14 Riemann invariant also
  setState("Down",4);

  // set shock slopes

  // (1): incident shock
  setShockSlope(1);

  // (2): reflected shock
  setShockSlope(2);

  // (3): mach stem
  setShockSlope(3);

  // set normal speed of the incident shock
  WWS = sqrt(pow((*WSh)(0,IP.at(0),ISH.at(0)),2) +
             pow((*WSh)(1,IP.at(0),ISH.at(0)),2));

  for(unsigned i=0; i<20; i++) { x.at(i) = xi.at(i); }

  computeTPstates.callCoUTP(xi,R14,R23,DXR14,DYR14,WWS);
  x = computeTPstates.getState();
  R14 = computeTPstates.getR();
  bool ifail = computeTPstates.getiFail();

  if(ifail) {
   xi.at(19)=-xi.at(19);
   computeTPstates.callCoUTP(xi,R14,R23,DXR14,DYR14,WWS);
   ifail = computeTPstates.getiFail();
   if (ifail) { 
    cout << "FixStateSps::error => triple point doesn't converge" << endl; }
   else {
    x = computeTPstates.getState();
    R14 = computeTPstates.getR();  }
  }
*/

  cout << "FixStateSps::error => TP (Triple point) scheme not";
  cout << " implemented yet\n";
  exit(1);
}

//--------------------------------------------------------------------------//

void FixStateSps::fixQPspecPoints(unsigned ISPPNTS)
{
  cout << "FixStateSps::error => QP (Quadruple point) scheme not";
  cout << " implemented yet\n";
  exit(1);
}

//--------------------------------------------------------------------------//

void FixStateSps::fixRRXspecPoints(unsigned ISPPNTS)
{
  cout << "FixStateSps::error => RRX (Regular reflection along x) scheme not";
  cout << " implemented yet\n";
  exit(1);
}

//--------------------------------------------------------------------------//

void FixStateSps::fixEPspecPoints(unsigned ISPPNTS)
{
  cout << "FixStateSps::error => EP (End point in a supersonic zone) scheme not";
  cout << " implemented yet\n";
  exit(1);
}

//--------------------------------------------------------------------------//

void FixStateSps::fixCspecPoints(unsigned ISPPNTS)
{
  cout << "FixStateSps::error => C (Connection between two shocks) scheme not";
  cout << " implemented yet\n";
  exit(1);
}

//--------------------------------------------------------------------------//

void FixStateSps::setState(string zone, unsigned nbState)
{
/*  long double help;
  Array3D <long double>* ZroeSh = 
   new Array3D <long double>((*ndof),(*nshmax),(*npshmax));

  if      (zone=="Up")   { ZroeSh = ZroeShu; }
  else if (zone=="Down") { ZroeSh = ZroeShd; }

  // set the start index to fill x vector which the state is referred to
  startIndex = (nbState-1) * 4;
  // set the shock indeces which the state is referred to
  ip = IP.at(nbState-1);
  ish = ISH.at(nbState-1);

  xi.at(startIndex+3) = (*ZroeSh)(3,ip,ish)/(*ZroeSh)(0,ip,ish); // y component
  xi.at(startIndex+2) = (*ZroeSh)(2,ip,ish)/(*ZroeSh)(0,ip,ish); // x component
  xi.at(startIndex+0) = (*ZroeSh)(1,ip,ish) * (*ZroeSh)(1,ip,ish); // density
  help = pow((*ZroeSh)(2,ip,ish),2) + pow((*ZroeSh)(3,ip,ish),2);
  xi.at(startIndex+1) = (*gm1)/(*gam) * ((*ZroeSh)(0,ip,ish) * 
                       (*ZroeSh)(1,ip,ish) - 0.5*help); // pressure

  // for the states (3) and (4) compute Riemann invariants 
  if ((nbState==3) || (nbState==4)) {
   dx = (*vShNor)(0,ip,ish);
   dy = (*vShNor)(0,ip,ish);
   if (nbState==3) { 
    R23 = sqrt((*gam)*xi.at(startIndex+1)/xi.at(startIndex+0)) +
          (*gm1) * 0.5 * (xi.at(startIndex+2) * dx + xi.at(startIndex+2) * dy);
   }
   else if (nbState==4) {
    R14 = sqrt((*gam)*xi.at(startIndex+1)/xi.at(startIndex+0)) +
          (*gm1) * 0.5 * (xi.at(startIndex+2) * dx + xi.at(startIndex+2) * dy);
   }
  }*/
}

//--------------------------------------------------------------------------//

void FixStateSps::setShockSlope(unsigned nbState)
{
/*  // set the start index to fill x vector which the state is referred to
  startIndex = (nbState-1) * 16;
  // set the shock indeces which the state is referred to
  ip = IP.at(nbState-1);
  ish = ISH.at(nbState-1);

  dx = (*vShNor)(0,ip,ish);
  dy = (*vShNor)(1,ip,ish);
  xi.at(startIndex) = atan2(-dx,dy);*/
}

//--------------------------------------------------------------------------//

void FixStateSps::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
{
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//--------------------------------------------------------------------------//

void FixStateSps::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDim();
  XYShu = new Array3D <double> (PhysicsInfo::getnbDim(),
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() + 
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim();
  XYShd = new Array3D <double> (PhysicsInfo::getnbDim(),
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZroeShu = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(), 
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void FixStateSps::freeArray()
{
  delete ZroeShu; delete ZroeShd;
  delete XYShu;   delete XYShd;
}

//--------------------------------------------------------------------------//

void FixStateSps::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void FixStateSps::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nSpecPoints =
      PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  nShockPoints =
      PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
      PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  typeSh =
      PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  typeSpecPoints =
      PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  XYSh = PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  SHinSPPs =
      PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
  ZroeShuOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZroeShdOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

