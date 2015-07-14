// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDps4Ar.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MoveDps4Ar, StateUpdater> moveDpsArProv("MoveDps4Ar");

//----------------------------------------------------------------------------//

MoveDps4Ar::MoveDps4Ar(const std::string& objectName) :
 StateUpdater(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDps4Ar::~MoveDps4Ar()
{
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::setup()
{
  LogToScreen(VERBOSE,"MoveDps4Ar::setup() => start\n");

  LogToScreen(VERBOSE,"MoveDps4Ar::setup() => end\n");
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::unsetup()
{
  LogToScreen(VERBOSE,"MoveDps4Ar::unsetup()\n");
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::update()
{
  LogToScreen(INFO,"MoveDps4Ar::update()\n");

  setMeshData();

  setPhysicsData();

  setAddress();

  double ZRHO;
  double RHOHF;

  logfile.Open(getClassName().c_str());

  WShMax = 0;

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
    RHOHF = RHOHF*ZRHO;
    RHOHF = 0;

    ro = ZRHO*ZRHO;
    help = pow((*ZroeSh)((*IX),IV,ISH),2) + pow((*ZroeSh)((*IY),IV,ISH),2);
    p = PhysicsInfo::getGm1()/PhysicsInfo::getGam() *
        (ZRHO * (*ZroeSh)((*IE),IV,ISH) - 0.5 * help * RHOHF);
    a = sqrt(PhysicsInfo::getGam()*p/ro);

    dum = 0;
    for(unsigned K=0; K<2; K++) { dum = dum + pow((*WSh)(K,IV,ISH),2); }
    WShMod = sqrt(dum);

    logfile("Shock point n. ", iShockPoint," speed: ", WShMod, "\n");

    if (WShMod>WShMax) { WShMax = WShMod; }

    dum = (MeshData::getInstance().getSHRELAX()) *
          (MeshData::getInstance().getDXCELL())  *
          (MeshData::getInstance().getSNDMIN()) /(a+WShMod);

    if(dt>dum) { dt = dum; }
   }
  }

  logfile("DT max: ", dt, "\n");

  // if the connectivity option is set to true and 
  // if the freezed connectivity ratio is less than the value
  // chosen in the input.case, the freezed connectivity bool variable
  // is set to true
  if(MeshData::getInstance().freezedConnectivityOption()) {

   // open file storing freezedConnectivityRatioValues
   ofstream memFCR("FreezedConnectivityRatioValues.dat",ios::app);
   memFCR << WShMax*dt/MeshData::getInstance().getSNDMIN() << endl;
   memFCR.close();

   if(WShMax*dt/MeshData::getInstance().getSNDMIN() <
      MeshData::getInstance().getFreezedAdimConnectivityRatio())
   {
    MeshData::getInstance().setFreezedConnectivity(true);
    logfile("\n (!) Freezed connectivity\n");
    logfile("Freezed connectivity ratio = ",
            WShMax*dt/MeshData::getInstance().getSNDMIN(),"\n\n");
    LogToScreen(INFO,"MoveDps4Ar::warning => freezed connectivity\n");
   }
  }

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

  // de-allocate ZRoeSh
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() + 
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeSh = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                 PhysicsInfo::getnbShPointsMax(),
                                 PhysicsInfo::getnbShMax(),
                                 &zroeVect->at(start));
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::freeArray()
{
  delete ZroeSh;
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void MoveDps4Ar::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  typeSh =
    PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  hf = PhysicsData::getInstance().getData <vector <double> > ("HF");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH"); 
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
