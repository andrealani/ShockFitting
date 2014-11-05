// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "WritingMeshSF/WriteSdwInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<WriteSdwInfo, WritingMesh> writeSdwInfoProv("WriteSdwInfo");

//--------------------------------------------------------------------------//

WriteSdwInfo::WriteSdwInfo(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteSdwInfo::~WriteSdwInfo()
{
  delete ZroeShu; delete ZroeShd;
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::setup()
{
  LogToScreen(VERBOSE,"WriteSdwInfo::setup() => start\n");

  LogToScreen(VERBOSE,"WriteSdwInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::unsetup()
{
  LogToScreen(VERBOSE,"WriteSdwInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::write()
{
  LogToScreen(INFO,"WriteSdwInfo::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  unsigned iShock;
  unsigned NSHE, IDUMMY;
  logfile.Open(getClassName().c_str());

  logfile("Opening file sh99.dat\n");

  file.open(string("sh99.dat").c_str());
  file.precision(15);

  file << (*nShocks) << "\n";
  logfile("nb. ",(*nShocks), " shock/discontinuities\n");

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   iShock = ISH+1;

   logfile("Shock/Discontinuity nb. ", iShock, "\n");
   file << nShockPoints->at(ISH) << " " << typeSh->at(ISH) << "\n";
   logfile("Kind of discontinuity: ", typeSh->at(ISH), "\n");
   logfile ("nb of points ",nShockPoints->at(ISH), "\n");

   for(unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for(unsigned IA=0; IA<(*ndim); IA++) { file << (*XYSh)(IA,I,ISH) << " "; }
    for(unsigned IA=0; IA<(*ndof); IA++) { file << (*ZroeShd)(IA,I,ISH) << " "; }
    for(unsigned IA=0; IA<(*ndof); IA++) { file << (*ZroeShu)(IA,I,ISH) << " "; }
    file << "\n";
   }
  }

  file << (*nSpecPoints) << "\n";
  logfile("\nnSpecPoints: ",(*nSpecPoints), "\n");

  for(unsigned ISPPNTS=0; ISPPNTS<(*nSpecPoints); ISPPNTS++) {
   file << typeSpecPoints->at(ISPPNTS) << "\n";
   logfile("Type Special Point: ",typeSpecPoints->at(ISPPNTS), "\n");

   // internal special point: triple point
   if ((*typeSpecPoints)[ISPPNTS]=="TP") {
    NSHE=4; IDUMMY=IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // internal special point: quad point
   else if ((*typeSpecPoints)[ISPPNTS]=="QP") {
    NSHE=5; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: regular reflection along x
   else if ((*typeSpecPoints)[ISPPNTS]=="RRX") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: wall point without reflection
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="WPNRX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: wall point without reflection
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="WPNRY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: inlet point
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="IPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: inlet point
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="IPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: outlet point
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="OPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: outlet point
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="OPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: end point in supersonic zone
   else if ((*typeSpecPoints)[ISPPNTS]=="EP") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: connection between two shocks
   else if ((*typeSpecPoints)[ISPPNTS]=="C") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    writeSHinSPPs(NSHE,ISPPNTS);
   }

   else { logfile("Condition not implemented");
          cout << "Condition not implemented\n";
          exit(1); }
  }

  file.close();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::writeSHinSPPs(unsigned NSHE, unsigned ISPPNTS)
{ 
  for (unsigned K=0; K<NSHE; K++) { 
   file << (*SHinSPPs)(0,K,ISPPNTS) << " "<< (*SHinSPPs)(1,K,ISPPNTS) << "\n";
   logfile((*SHinSPPs)(0,K,ISPPNTS)," ",(*SHinSPPs)(1,K,ISPPNTS),"\n" );
  }
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::setAddress()
{
  unsigned start;
  start = npoin->at(0)*(*ndof);
  ZroeShu = 
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
  start = npoin->at(0) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  ZroeShd = 
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void WriteSdwInfo::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  typeSpecPoints =
     PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  typeSh = PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
  SHinSPPs =
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

