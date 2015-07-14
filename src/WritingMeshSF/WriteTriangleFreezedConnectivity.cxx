// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteTriangleFreezedConnectivity.hh"
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
ObjectProvider<WriteTriangleFreezedConnectivity, WritingMesh> 
 writeTriangleFreezedConnectivityProv("WriteTriangleFreezedConnectivity");

//--------------------------------------------------------------------------//

WriteTriangleFreezedConnectivity::WriteTriangleFreezedConnectivity
(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteTriangleFreezedConnectivity::~WriteTriangleFreezedConnectivity()
{
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setup()
{
  LogToScreen(VERBOSE,"WriteTriangleFreezedConnectivity::setup() => start\n");

  LogToScreen(VERBOSE,"WriteTriangleFreezedConnectivity::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::unsetup()
{
  LogToScreen(VERBOSE,"WriteTriangleFreezedConnectivity::unsetup()\n");
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::write()
{
  LogToScreen(INFO,"WriteTriangleFreezedConnectivity:write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  // copy and rename the triangle files
  copyAnDrenameTriangleFiles();

  ilist = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                             PhysicsInfo::getnbShPointsMax();

  // M02M1 and M12M0 are filled with indeces that start
  // from 1 to NPOIN+2*NSHMAX*NPSHMAX+1
  M02M1->resize(ilist+1); // c++ indeces start from 0
  M12M0->resize(ilist+1); // c++ indeces start from 0

  TNPOIN = 0;

  // set map vector for nodcod
  setMapVectorForNodcod();

  icount = npoin->at(0);

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  ilist = TNPOIN;

  npoin->at(1) = ilist;

  // fill in the .node file
  writeFileNode();

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::copyAnDrenameTriangleFiles()
{
  string command;
  string command_comm;
  string fnameOld;
  stringstream backdir;

  fnameOld = fname->str();

  fname->str(string());

  unsigned nbDig = 0;
  unsigned dummyIstep = MeshData::getInstance().getIstep();

  while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
  *fname << setw(7-nbDig) << setfill('0') << left << string("na").c_str();
  *fname << MeshData::getInstance().getIstep();

  nbDig = 0;
  dummyIstep = MeshData::getInstance().getIstep()-1;
  while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
  backdir << setw(9-nbDig) << setfill('0') << left 
          << string("step").c_str() << MeshData::getInstance().getIstep()-1;

  // if the solution has been saved in the previous step
  // take the files inside the corresponding folder
  if(((MeshData::getInstance().getIstep()-1)-1)% MeshData::getInstance().getnbIbak()==0) {
   command_comm = "cp " + MeshData::getInstance().getResultsDir() + "/" + backdir.str() + "/" +fnameOld;
   command = command_comm + ".1.poly ./" + fname->str() + ".1.poly";
   system(command.c_str());
   command = command_comm + ".1.ele ./" + fname->str() + ".1.ele";
   system(command.c_str());
   command = command_comm + ".1.neigh ./" + fname->str() + ".1.neigh";
   system(command.c_str());
   command = command_comm + ".1.edge ./" + fname->str() + ".1.edge";
   system(command.c_str());
  }

  // else copy the triangle files stored in the execution directory 
  else if ((((MeshData::getInstance().getIstep()-1)-1)% MeshData::getInstance().getnbIbak()!=0) && (MeshData::getInstance().getIstep()-1)!=1) { 
   command = "cp " + fnameOld + ".1.poly " + fname->str() + ".1.poly";
   system(command.c_str());
   command = "cp " + fnameOld + ".1.ele " + fname->str() + ".1.ele";
   system(command.c_str());
   command = "cp " + fnameOld + ".1.neigh " + fname->str() + ".1.neigh";
   system(command.c_str());
   command = "cp " + fnameOld + ".1.edge " + fname->str() + ".1.edge"; 
   system(command.c_str());
  }

  else if ((MeshData::getInstance().getIstep()-1)==1) {
   command_comm = "cp " + MeshData::getInstance().getResultsDir() + "/step00001/";
   command = command_comm + "na00001.1.poly ./na00002.1.poly";
   system(command.c_str());
   command = command_comm + "na00001.1.ele ./na00002.1.ele";
   system(command.c_str());
   command = command_comm + "na00001.1.neigh ./na00002.1.neigh";
   system(command.c_str());
   command = command_comm + "na00001.1.edge ./na00002.1.edge";
   system(command.c_str());
  }

}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::writeFileNode()
{
  string dummyfile;

  // open node file
  dummyfile = fname->str()+".1.node";

  file = fopen(dummyfile.c_str(),"w");

  fprintf(file,"%u %s %u %s %u %s",ilist," ",
          PhysicsInfo::getnbDim()," ",*ndof," 1\n");

  // write mesh variables
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if(nodcod->at(IPOIN)>=0) {
    fprintf(file,"%u %s",M02M1->at(IPOIN+1),"  "); // c++ indeces start from 0
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file,"%.16f %s",(*XY)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%u %s",nodcod->at(IPOIN),"\n");
   }
  }

  icount = npoin->at(0);

  // write upstream status
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
    ++icount;
    if((*NodCodSh)(I,ISH)==10) {
     fprintf(file,"%u %s",M02M1->at(icount),"  ");
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file,"%0.16f %s",(*XYShu)(IA,I,ISH),"  ");}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%0.16f %s",(*ZRoeShu)(IA,I,ISH)," ");}
     fprintf(file,"%u %s",(*NodCodSh)(I,ISH),"\n");
    }
   }
  }

  // write downstream status
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
     icount++;
    if((*NodCodSh)(I,ISH)==10) {
     fprintf(file,"%u %s",M02M1->at(icount),"  ");
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file, "%0.16f %s",(*XYShd)(IA,I,ISH),"  ");}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%0.16f %s",(*ZRoeShd)(IA,I,ISH),"  ");}
     fprintf(file,"%u %s",(*NodCodSh)(I,ISH),"\n");
    }
   }
  }

  fclose(file);
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setMapVectorForNodcod()
{
  for (unsigned IPOIN=0; IPOIN< npoin->at(0); IPOIN++) {
   if (nodcod->at(IPOIN)>=0) { ++TNPOIN;
                               M02M1->at(IPOIN+1) = TNPOIN;
                               M12M0->at(TNPOIN) = IPOIN+1; }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setMapVectorForNodcodSh()
{
  for (unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
    ++icount;

    if ((*NodCodSh)(I,ISH)==10) { ++TNPOIN;
                                    M02M1->at(icount) = TNPOIN;
                                    M12M0->at(TNPOIN) = icount;}
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setAddress()
{
  unsigned start; unsigned totsize;
  start = 0;
  XY = new Array2D <double> (PhysicsInfo::getnbDim(), npoin->at(0),
                             &coorVect->at(start));
  Zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),npoin->at(0),
                               &zroeVect->at(start));
  totsize = nbfac->at(0) + 2 *
                           PhysicsInfo::getnbShMax() *
                           PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  celnod = new Array2D<int> ((*nvt), nelem->at(0), &celnodVect->at(start));
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZRoeShu = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() + 
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() ;
  XYShu = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim()  +
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim() ;
  XYShd = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::freeArray()
{
  delete XY; delete Zroe;
  delete bndfac; delete celnod;
  delete ZRoeShu; delete ZRoeShd; delete XYShu; delete XYShd; 
  delete NodCodSh;
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem =MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbfacSh = MeshData::getInstance().getData<unsigned>("NBFACSH");
  caddholes =
    MeshData::getInstance().getData <vector<double> > ("CADDholes");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  fname = MeshData::getInstance().getData <stringstream> ("FNAME");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
}

//--------------------------------------------------------------------------//

void WriteTriangleFreezedConnectivity::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
