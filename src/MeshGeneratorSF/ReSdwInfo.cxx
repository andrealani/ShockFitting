// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReSdwInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
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
}

//--------------------------------------------------------------------------//

void ReSdwInfo::unsetup()
{
  LogToScreen(VERBOSE, "ReSdwInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReSdwInfo::generate()
{
  LogToScreen(INFO, "ReadSdwInfo::generate()\n");

  logfile.Open(getClassName());
  // check if the sh00 file is already opened
  if(file.is_open()) { 
   cout << "ReSdwInfo::warning => file sh00.dat seems to be already opened\n";}

  file.open(getInputFiles().c_str());

  setMeshData();
  setPhysicsData();

  setAddress();
  setSize();

  readShockInfo(); 

  file.close();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReSdwInfo::generate(string processingFile)
{
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
  for (unsigned ISH=0; ISH < PhysicsInfo::getnbShMax(); ISH++) {
    for (unsigned K=0; K < PhysicsInfo::getnbShPointsMax(); K++)
     {(*NodCodSh)(K,ISH) = -99;}
  }
  logfile("Open sh00.dat\n");

  file >> (*nShocks);
  logfile("Found n. ",(*nShocks)," shock/discontinuities\n");

  file.precision(18);

  for (unsigned ISH=0; ISH < (*nShocks); ISH++) {
   unsigned iShock = ISH+1; // c++ indeces start from 0
   logfile("Shock/Discontinuity n. ",iShock,"\n");
   file >> (*nShockPoints)[ISH] >> (*typeSh)[ISH] ;
   logfile("Kind of discontinuity: ",(*typeSh)[ISH],"\n");
   logfile("n. of points ",(*nShockPoints)[ISH],"\n");

   (*nShockEdges)[ISH] = (*nShockPoints)[ISH]-1;
   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I <PhysicsInfo::getnbDim(); I++) {
     file >> (*XYSh)(I,K,ISH);
     logfile((*XYSh)(I,K,ISH)," ");
    }
    for (unsigned I=0; I < (*ndof); I++) {
     file >> (*ZRoeShd)(I,K,ISH);
     logfile((*ZRoeShd)(I,K,ISH)," ");
    }
    for (unsigned I=0; I < (*ndof); I++) {
     file >> (*ZRoeShu)(I,K,ISH);
     logfile((*ZRoeShu)(I,K,ISH)," ");
    }
    logfile("\n");
    (*NodCodSh)(K,ISH) = 10;
   }

   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndof); I++) {
     (*ZRoeShuOld)(I,K,ISH) = (*ZRoeShu)(I,K,ISH);
     (*ZRoeShdOld)(I,K,ISH) = (*ZRoeShd)(I,K,ISH);
    }
   }
  }

  IDUMMY = 0;

  file >> (*nSpecPoints);
  typeSpecPoints->resize((*nSpecPoints));
  logfile("nSpecPoints: ",*nSpecPoints,"\n");
  for (unsigned ISPPNTS=0; ISPPNTS < (*nSpecPoints); ISPPNTS++) {
   file >> (*typeSpecPoints)[ISPPNTS];
   logfile("Type Special Point: ",(*typeSpecPoints)[ISPPNTS],"\n");

   // internal special point: triple point
   if ((*typeSpecPoints)[ISPPNTS]=="TP") {
    NSHE=4; IDUMMY=IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // internal special point: quad point
   else if ((*typeSpecPoints)[ISPPNTS]=="QP") {
    NSHE=5; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: regular reflection along x
   else if ((*typeSpecPoints)[ISPPNTS]=="RRX") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: wall point without reflection
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="WPNRX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: wall point without reflection
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="WPNRY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: inlet point
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="IPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: inlet point
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="IPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: outlet point
   // floating along X direction
   else if ((*typeSpecPoints)[ISPPNTS]=="OPX") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: outlet point
   // floating along Y direction
   else if ((*typeSpecPoints)[ISPPNTS]=="OPY") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }   

   // boundary special point: end point in supersonic zone
   else if ((*typeSpecPoints)[ISPPNTS]=="EP") {
    NSHE=1; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   // boundary special point: connection between two shocks
   else if ((*typeSpecPoints)[ISPPNTS]=="C") {
    NSHE=2; IDUMMY = IDUMMY+NSHE;
    setSHinSPPs(NSHE,ISPPNTS);
   }

   else { logfile("Condition not implemented");
          cout << "Condition not implemented\n";
          exit(1); }
  }
 
  // check conditions on special points
  if (IDUMMY != 2 * (*nShocks)) {
   logfile("Nof conditions on special points not correct");
   unsigned ask = 2 * (*nShocks);
   logfile("Nof asked conditions: ",ask,"\n");
   logfile("Nof forced conditions: ",IDUMMY,"\n");
   cout << "Nof conditions on special points not correct\n";
   exit(1);
  }
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setSHinSPPs(unsigned NSHE, unsigned ISPPNTS)
{
  for (unsigned K=0; K<NSHE; K++) {
   file >> (*SHinSPPs)(0,K,ISPPNTS) >> (*SHinSPPs)(1,K,ISPPNTS);
   logfile((*SHinSPPs)(0,K,ISPPNTS)," ",(*SHinSPPs)(1,K,ISPPNTS),"\n" );
  }
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setSize()
{
  nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  typeSpecPoints->resize((*nSpecPoints));
  typeSh->resize(PhysicsInfo::getnbShMax());
  XYSh->resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax());
  ZRoeShuOld->resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  ZRoeShdOld->resize((*ndof),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  SHinSPPs->resize(2,5,PhysicsInfo::getnbSpecPointsMax());
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setAddress()
{
  unsigned start;
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0)*PhysicsInfo::getnbDofMax();
  ZRoeShu = 
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroe->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() + 
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd = 
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroe->at(start));
}

//--------------------------------------------------------------------------//

void ReSdwInfo::freeArray()
{
  delete NodCodSh; delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nodcod = MeshData::getInstance().getData <vector <int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void ReSdwInfo::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
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
  ZRoeShuOld = 
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZRoeShdOld = 
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
  SHinSPPs = 
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
