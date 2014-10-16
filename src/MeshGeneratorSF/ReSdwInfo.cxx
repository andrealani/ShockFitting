// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReSdwInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "MathTools/Array3D.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReSdwInfo, MeshGenerator> readSdwInfoProv("ReSdwInfo");

//--------------------------------------------------------------------------//

ReSdwInfo::ReSdwInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

ReSdwInfo::~ReSdwInfo()
{
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setup()
{
  LogToScreen(VERBOSE, "ReSdwInfo::setup() => start\n");

  setMeshData();
  setPhysicsData();

  logfile.Open(getClassName());
}

//--------------------------------------------------------------------------//

void ReSdwInfo::unsetup()
{
  LogToScreen(VERBOSE, "ReSdwInfo::unsetup()\n");

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReSdwInfo::generate()
{
  LogToScreen(INFO, "ReadSdwInfo::generate()\n");

  file.open(getInputFiles().c_str());

  setAddress();
  setSize();

  readShockInfo(); 

  file.close();
}

//--------------------------------------------------------------------------//

std::string ReSdwInfo::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

void ReSdwInfo::readShockInfo()
{
  unsigned NSHE, IDUMMY;

  // initialize NODCODSH which is part of NODCOD
  // If the code -99 is used this means no shock point
  // If the code 10 is used this means shock point
  for (unsigned ISH=0; ISH < (*nshmax); ISH++) {
    for (unsigned K=0; K < (*npshmax); K++) {(*r_NodCodSh)(K,ISH) = -99;}
  }
  logfile("Open sh00.dat\n");

  file >> (*r_nShocks);
  logfile("Found n. ",(*r_nShocks)," shock/discontinuities\n");

  for (unsigned ISH=0; ISH < (*r_nShocks); ISH++) {
   unsigned iShock = ISH+1; // c++ indeces start from 0
   logfile("Shock/Discontinuity n. ",iShock,"\n");
   file >> (*r_nShockPoints)[ISH] >> (*r_typeSh)[ISH] ;
   logfile("Kind of discontinuity: ",(*r_typeSh)[ISH],"\n");
   logfile("n. of points ",(*r_nShockPoints)[ISH],"\n");

   (*r_nShockEdges)[ISH] = (*r_nShockPoints)[ISH]-1;
   for (unsigned K=0; K < (*r_nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndim); I++) {
     file >> (*r_XYSh)(I,K,ISH);
     logfile((*r_XYSh)(I,K,ISH)," ");
    }
    for (unsigned I=0; I < (*ndof); I++) {
     file >> (*r_ZRoeShd)(I,K,ISH);
     logfile((*r_ZRoeShd)(I,K,ISH)," ");
    }
    for (unsigned I=0; I < (*ndof); I++) {
     file >> (*r_ZRoeShu)(I,K,ISH);
     logfile((*r_ZRoeShu)(I,K,ISH)," ");
    }
    logfile("\n");
    (*r_NodCodSh)(K,ISH) = 10;
   }

   for (unsigned K=0; K < (*r_nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndof); I++) {
     (*r_ZRoeShuOld)(I,K,ISH) = (*r_ZRoeShu)(I,K,ISH);
     (*r_ZRoeShdOld)(I,K,ISH) = (*r_ZRoeShd)(I,K,ISH);
    }
   }
  }

  IDUMMY = 0;

  file >> (*r_nSpecPoints);
  r_typeSpecPoints->resize((*r_nSpecPoints));
  logfile("nSpecPoints: ",*r_nSpecPoints,"\n");
  for (unsigned ISPPNTS=0; ISPPNTS < (*r_nSpecPoints); ISPPNTS++) {
   file >> (*r_typeSpecPoints)[ISPPNTS];
   logfile("Type Special Point: ",(*r_typeSpecPoints)[ISPPNTS],"\n");

   // internal special point: triple point
   if ((*r_typeSpecPoints)[ISPPNTS]=="TP") {
    NSHE=4; IDUMMY=IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // internal special point: quad point
   else if ((*r_typeSpecPoints)[ISPPNTS]=="QP") {
    NSHE=5; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: regular reflection along x
   else if ((*r_typeSpecPoints)[ISPPNTS]=="RRX") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: wall point without reflection
   // floating along X direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="WPNRX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: wall point without reflection
   // floating along Y direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="WPNRY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: inlet point
   // floating along X direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="IPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: inlet point
   // floating along Y direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="IPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: outlet point
   // floating along X direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="OPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    } 
   }

   // boundary special point: outlet point
   // floating along Y direction
   else if ((*r_typeSpecPoints)[ISPPNTS]=="OPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }   

   // boundary special point: end point
   else if ((*r_typeSpecPoints)[ISPPNTS]=="EP") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   // boundary special point: connection between two shocks
   else if ((*r_typeSpecPoints)[ISPPNTS]=="C") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    for (unsigned K=0; K<NSHE; K++) {
     file >> (*r_SHinSPPs)(0,K,ISPPNTS) >> (*r_SHinSPPs)(1,K,ISPPNTS);
     logfile((*r_SHinSPPs)(0,K,ISPPNTS),(*r_SHinSPPs)(1,K,ISPPNTS),"\n" );
    }
   }

   else {
    logfile("Condition not implemented");
    cout << "Condition not implemented\n";
    exit(1);
   }
  }
 
  // check conditions on special points
  if (IDUMMY != 2 * (*r_nShocks)) {
   logfile("Nof conditions on special points not correct");
   unsigned ask = 2 * (*r_nShocks);
   logfile("Nof asked conditions: ",ask,"\n");
   logfile("Nof forced conditions: ",IDUMMY,"\n");
   cout << "Nof conditions on special points not correct\n";
   exit(1);
  }
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nodcod = MeshData::getInstance().getData <vector <int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nspmax = PhysicsData::getInstance().getData <unsigned> ("NSPMAX");
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");

  r_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  r_nShockPoints = 
       PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  r_nShockEdges = 
       PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  r_nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  r_typeSpecPoints = 
       PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  r_typeSh = PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  r_XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
  r_ZRoeShuOld = 
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  r_ZRoeShdOld = 
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
  r_SHinSPPs = 
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setAddress()
{
  unsigned start;
  start = (*npoin);
  r_NodCodSh = new Array2D <int> ((*npshmax),(*nshmax),&nodcod->at(start));
  start = (*npoin)*(*ndof);
  r_ZRoeShu = new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
  start = (*npoin) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  r_ZRoeShd = new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroe->at(start));
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setSize()
{
  r_nShockPoints->resize((*nshmax));
  r_nShockEdges->resize((*nshmax));
  r_typeSpecPoints->resize((*r_nSpecPoints));
  r_typeSh->resize((*nshmax));
  r_XYSh->resize((*ndim),(*npshmax),(*nshmax));
  r_ZRoeShuOld->resize((*ndofmax),(*npshmax),(*nshmax));
  r_ZRoeShdOld->resize((*ndof),(*npshmax),(*nshmax));
  r_SHinSPPs->resize(2,5,(*nspmax));
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
