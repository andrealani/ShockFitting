// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/RdDps.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<RdDps, Remeshing> redistributeShPointsProv("RdDps");

//--------------------------------------------------------------------------//

RdDps::RdDps(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

RdDps::~RdDps()
{
  delete ZroeShu; delete ZroeShd;
}

//--------------------------------------------------------------------------//

void RdDps::setup()
{
  LogToScreen(VERBOSE,"RdDps::setup() => start\n");

  LogToScreen(VERBOSE,"RdDps::setup() => end\n");
}

//--------------------------------------------------------------------------//

void RdDps::unsetup()
{
  LogToScreen(VERBOSE,"RdDps::unsetup()\n");
}

//--------------------------------------------------------------------------//

void RdDps::remesh()
{
  LogToScreen(INFO,"RdDps::remesh()\n");

  setMeshData();
  setPhysicsData();

  setAddress();
/*
ifstream var;
stringstream pathvar;
pathvar.str(string());
if(MeshData::getInstance().getIstep()<10){
//pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.0/tests/CircularCylinder-Unibas_inv_M20_coarse_N/step0000"<<MeshData::getInstance().getIstep()<<"/Var/rddps.var";
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step0000"<<MeshData::getInstance().getIstep()<<"/Var/rddps.var";
}
else if (MeshData::getInstance().getIstep()>=10 &&
         MeshData::getInstance().getIstep()<100){
//pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.0/tests/CircularCylinder-Unibas_inv_M20_coarse_N/step000"<<MeshData::getInstance().getIstep()<<"/Var/rddps.var";
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step000"<<MeshData::getInstance().getIstep()<<"/Var/rddps.var";
}

else if (MeshData::getInstance().getIstep()>=100 &&
         MeshData::getInstance().getIstep()<1000){
pathvar << "/students/st_13_14/deamicis/nobackup/UnDiFi-2D-v2.1/tests/CircularCylinder_VKI_inv_N-N2_E2_LRD/step00"<<MeshData::getInstance().getIstep()<<"/Var/rddps.var";
}


string path = pathvar.str();
var.open(path.c_str());

if(var.fail()) { cout << "Step000" << MeshData::getInstance().getIstep() << "Failed opening rddps.var" << endl;
}

  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for(unsigned k=0;k<(*ndof);k++) { var >> (*ZroeShu)(k,I,ISH);}
    for(unsigned k=0;k<(*ndof);k++) { var >> (*ZroeShd)(k,I,ISH);}
    for(unsigned k=0;k<2;k++) { var >> (*XYSh)(k,I,ISH);}
}}
var.close();
*/

  unsigned iShPoint;
  double dum;
  unsigned ileMin, ileMax;
  double lenRelMin, lenRelMax;

  logfile.Open(getClassName().c_str());

  ShEdgeLgth.resize(PhysicsInfo::getnbShPointsMax());

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   ileMin = 0; ileMax = 0;
   lenRelMin = 1; lenRelMax = 1;

   // compute length of shock edge
   for(unsigned IV=0; IV<nShockEdges->at(ISH); IV++) {
    ShEdgeLgth.at(IV) = 0;
    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     ShEdgeLgth.at(IV) = ShEdgeLgth.at(IV) + 
                         pow(((*XYSh)(I,IV,ISH)-(*XYSh)(I,IV+1,ISH)),2);
    }
    ShEdgeLgth.at(IV) = sqrt(ShEdgeLgth.at(IV));
    dum = ShEdgeLgth.at(IV)/(MeshData::getInstance().getDXCELL());

    if (dum<0.5) {
     if (lenRelMin>dum) { lenRelMin = dum;
                          ileMin = IV+1;   }
    }
    if (dum>1.5) {
     if (lenRelMax<dum) { lenRelMax = dum;
                          ileMax = IV+1; } 
    }
   }

   if (nShockPoints->at(ISH)<=3) { ileMin = 0;
                                   ileMax = 0; }


   if (ileMin != 0) {

    cout << "\n\n\nRdDps::warning => ileMin!=0 this part has never been tested\n";
    cout << "                  if something goes wrong please check the IV index\n";
    cout << "                  and the output file ileMinLog.check\n\n\n";

    logfile("Before\n ");
    for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
     iShPoint = IV+1;
     logfile(iShPoint," ", (*ZroeShd)(0,IV,ISH)," ",(*ZroeShd)(1,IV,ISH),"\n");
    }
    unsigned npc = ileMin;

    if (ShEdgeLgth.at(ileMin-2)>ShEdgeLgth.at(ileMin)  )   { npc = ileMin+1; }
    if (ileMin==1)                                         { npc = 2; }
    if (ileMin==nShockEdges->at(ISH))                      { npc = ileMin; }

    for(unsigned IV=npc; IV<nShockPoints->at(ISH); IV++) {
     for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
      (*XYSh)(I,IV-1,ISH) = (*XYSh)(I,IV,ISH); // c++ indeces start from 0
     }
     for(unsigned I=0; I<(*ndof); I++) {
      (*ZroeShu)(I,IV-1,ISH) = (*ZroeShu)(I,IV,ISH); // c++ indeces start from 0
      (*ZroeShd)(I,IV-1,ISH) = (*ZroeShd)(I,IV,ISH); // c++ indeces start from 0
     }
    }

    nShockPoints->at(ISH) = nShockPoints->at(ISH)- 1;
    nShockEdges->at(ISH) = nShockEdges->at(ISH)- 1;

    logfile("After\n");
    for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
     iShPoint = IV+1;
     logfile(iShPoint," ", (*ZroeShd)(0,IV,ISH)," ",(*ZroeShd)(1,IV,ISH),"\n");
    }

    FILE* ileMinLog;
    ileMinLog = fopen("IleMinLog.check","w");

     for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
      for(unsigned k=0;k<(*ndof);k++) {
       fprintf(ileMinLog,"%32.16F %s",(*ZroeShd)(k,IV,ISH)," ");}
      for(unsigned k=0;k<(*ndof);k++) {
       fprintf(ileMinLog,"%32.16F %s",(*ZroeShu)(k,IV,ISH)," ");}
        for(unsigned k=0;k<2;k++) { 
       fprintf(ileMinLog,"%32.16F %s",(*XYSh)(k,IV,ISH)," ");}
     }

    fclose(ileMinLog);

   }

   if (ileMax !=0) { 
    logfile("Before\n");
    for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
     iShPoint = IV+1;
     logfile(iShPoint," ", (*ZroeShd)(0,IV,ISH)," ",(*ZroeShd)(1,IV,ISH),"\n");
    }

    unsigned npi = ileMax;
    for(unsigned IV=nShockPoints->at(ISH)-1;IV>=npi;IV--) {
     for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
      (*XYSh)(I,IV+1,ISH) = (*XYSh)(I,IV,ISH); // c++ indeces start from 0
     }
     for(unsigned I=0; I<(*ndof); I++) {
      (*ZroeShu)(I,IV+1,ISH) = (*ZroeShu)(I,IV,ISH); // c++ indeces start from 0
      (*ZroeShd)(I,IV+1,ISH) = (*ZroeShd)(I,IV,ISH); // c++ indeces start from 0
     }
    }
    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     // c++ indeces start from 0 
     (*XYSh)(I,npi,ISH) = 0.5*((*XYSh)(I,npi-1,ISH)+(*XYSh)(I,npi+1,ISH));
    }
    for(unsigned I=0; I<(*ndof); I++) {
     // c++ indeces start from 0
     (*ZroeShu)(I,npi,ISH) = 
        0.5 * ( (*ZroeShu)(I,npi-1,ISH) + (*ZroeShu)(I,npi+1,ISH) );
     (*ZroeShd)(I,npi,ISH) = 
        0.5 * ( (*ZroeShd)(I,npi-1,ISH) + (*ZroeShd)(I,npi+1,ISH) );
    }

    nShockPoints->at(ISH) = nShockPoints->at(ISH)+1;
    nShockEdges->at(ISH) = nShockEdges->at(ISH)+1;

    logfile("After\n");
    for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
     iShPoint = IV+1;
     logfile(iShPoint," ", (*ZroeShd)(0,IV,ISH)," ",(*ZroeShd)(1,IV,ISH),"\n");
    }
   }
  }

  logfile.Close();

FILE* output;
output =fopen("CheckC/rddps.check","w");

for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
    for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
     for(unsigned k=0;k<(*ndof);k++) {
      fprintf(output,"%32.16F %s",(*ZroeShd)(k,IV,ISH)," ");}
     for(unsigned k=0;k<(*ndof);k++) {
      fprintf(output,"%32.16F %s",(*ZroeShu)(k,IV,ISH)," ");}
         for(unsigned k=0;k<2;k++) {
      fprintf(output,"%32.16F %s",(*XYSh)(k,IV,ISH)," ");}
}}

fclose(output);
}

//--------------------------------------------------------------------------//

void RdDps::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZroeShu = new Array3D <double> 
              (PhysicsInfo::getnbDofMax(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax(),
               &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeShd = new Array3D <double> 
              (PhysicsInfo::getnbDofMax(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax(),
               &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void RdDps::setMeshData()
{
  zroeVect = MeshData::getInstance().getData <vector <double> >("ZROE");
  npoin  = MeshData::getInstance().getData <vector<unsigned> >("NPOIN"); 
}

//--------------------------------------------------------------------------//

void RdDps::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = 
    PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges = 
    PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  XYSh = 
    PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
