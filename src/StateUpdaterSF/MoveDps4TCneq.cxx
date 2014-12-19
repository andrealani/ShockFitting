// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDps4TCneq.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MoveDps4TCneq, MoveDps> moveDpsTCneqProv("MoveDps4TCneq");

//----------------------------------------------------------------------------//

MoveDps4TCneq::MoveDps4TCneq(const std::string& objectName) :
 MoveDps(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDps4TCneq::~MoveDps4TCneq()
{
}

//----------------------------------------------------------------------------//

void MoveDps4TCneq::setup()
{
  LogToScreen(VERBOSE,"MoveDps4TCneq::setup() => start\n");

  LogToScreen(VERBOSE,"MoveDps4TCneq::setup() => end\n");
}

//----------------------------------------------------------------------------//

void MoveDps4TCneq::unsetup()
{
  LogToScreen(VERBOSE,"MoveDps4TCneq::unsetup()\n");
}

//----------------------------------------------------------------------------//

void MoveDps4TCneq::update()
{
  LogToScreen(INFO,"MoveDps4TCneq::update()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

/*
ifstream var;
stringstream pathvar;
pathvar.str(string());
if(MeshData::getInstance().getIstep()<10){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step0000"<<MeshData::getInstance().getIstep()<<"/Var/mvdps.var";
}
else if (MeshData::getInstance().getIstep()>=10 &&
         MeshData::getInstance().getIstep()<100){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step000"<<MeshData::getInstance().getIstep()<<"/Var/mvdps.var";
}
else if (MeshData::getInstance().getIstep()>=100 &&
         MeshData::getInstance().getIstep()<1000){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step00"<<MeshData::getInstance().getIstep()<<"/Var/mvdps.var";
}


string path = pathvar.str();
var.open(path.c_str());

if(var.fail()) { cout << "Step000" << MeshData::getInstance().getIstep() << "Failed opening mvdps.var" << endl;
}

  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for(unsigned k=0;k<(*ndof);k++) { var >> (*ZroeSh)(k,I,ISH);}
    for(unsigned k=0;k<2;k++) { var >> (*XYSh)(k,I,ISH);}
    for(unsigned k=0;k<2;k++) { var >> (*WSh)(k,I,ISH);}
}}
var.close();


*/

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
    help = pow((*ZroeSh)((*IX),IV,ISH),2)+pow((*ZroeSh)((*IY),IV,ISH),2);
    p = ((*gref)-1.0)/(*gref) * (ZRHO * (*ZroeSh)((*IE),IV,ISH) - 
         0.50 * help - RHOHF - (*ZroeSh)((*IEV),IV,ISH) * ZRHO);
    a = sqrt((*gref)*p/ro);

    dum = 0;
    for(unsigned K=0; K<2; K++) { dum = dum + pow((*WSh)(K,IV,ISH),2); }
    WShMod = sqrt(dum);

    logfile("Shock point n. ", iShockPoint," speed: ", WShMod, "\n");

    dum = (MeshData::getInstance().getSHRELAX()) *
          (MeshData::getInstance().getDXCELL())  *
          (MeshData::getInstance().getSNDMIN()) /(a+WShMod);
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

FILE* output;
output =fopen("CheckC/mvdps.check","w");

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    for(unsigned I=0; I<2; I++) {
     fprintf(output,"%32.16F %s",(*XYSh)(I,IV,ISH), " ");}
     fprintf(output,"%s","\n");}}
fclose(output);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
