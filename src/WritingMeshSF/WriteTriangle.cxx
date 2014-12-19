// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteTriangle.hh"
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
ObjectProvider<WriteTriangle, WritingMesh> writeTriangleProv("WriteTriangle");

//--------------------------------------------------------------------------//

WriteTriangle::WriteTriangle(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteTriangle::~WriteTriangle()
{
}

//--------------------------------------------------------------------------//

void WriteTriangle::setup()
{
  LogToScreen(VERBOSE, "WriteTriangle::setup() => start\n");

  LogToScreen(VERBOSE, "WriteTriangle::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteTriangle::unsetup()
{
  LogToScreen(VERBOSE, "WriteTriangle::unsetup()\n");
}

//--------------------------------------------------------------------------//

void WriteTriangle::write()
{
  LogToScreen(INFO, "WriteTriangle::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  string dummyfile;

  fname->str(string());

  unsigned nbDig = 0;
  unsigned dummyIstep = MeshData::getInstance().getIstep();

  while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
  *fname << setw(7-nbDig) << setfill('0') << left << string("na").c_str();
  *fname << MeshData::getInstance().getIstep();

  // write node file
  dummyfile = fname->str()+".node";

  file = fopen(dummyfile.c_str(),"w");

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

  fprintf(file,"%u %s %u %s %u %s",ilist," ",PhysicsInfo::getnbDim()," ",*ndof," 1\n");
  
  // write mesh points coordinates and states on Triangle file
  writeMeshVariables();

  icount = npoin->at(0);

  // write upstream shock points coordinates and states on Triangle file
  writeUpstreamStatus();

  // write downstream shock points coordinates and states on Triangle file
  writeDownstreamStatus();

  fclose(file);

  // write poly file
  dummyfile = fname->str()+".poly";
  file = fopen(dummyfile.c_str(),"w");

  fprintf(file,"%s %u %s", "0 ",PhysicsInfo::getnbDim()," 0 1 \n");
  
  // set map vector for bndfac and write bndfac on poly file
  writeBndfac();

  // compute number of holes
  computenbHoles();

  fclose(file);
}

//--------------------------------------------------------------------------//

void WriteTriangle::setMapVectorForNodcod()
{
  for (unsigned IPOIN=0; IPOIN< npoin->at(0); IPOIN++) {
   if (nodcod->at(IPOIN)>=0) { ++TNPOIN;
                               M02M1->at(IPOIN+1) = TNPOIN;
                               M12M0->at(TNPOIN) = IPOIN+1; } 
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::setMapVectorForNodcodSh()
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

void WriteTriangle::writeMeshVariables()
{
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
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeUpstreamStatus()
{
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
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeDownstreamStatus()
{
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
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeBndfac()
{
  int IBC;

  ICHECK = 0;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(2,IFACE);
   if(IBC>0) { ++ICHECK; }
  }
  fprintf(file,"%u %s",ICHECK,"  1\n");
  ICHECK=0;

  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) { 
   IBC = (*bndfac)(2,IFACE);

   if(IBC>0) {
    ++ICHECK; 
    fprintf(file,"%u %s",ICHECK,"  ");
    for(unsigned IA=0; IA<2; IA++) {
     fprintf(file,"%u %s",M02M1->at((*bndfac)(IA,IFACE)),"  ");
    }
   fprintf(file,"%i %s",IBC,"\n");
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::computenbHoles()
{
  nHoles = 0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   nHoles = nHoles + nShockPoints->at(ISH)-2;
  }

  fprintf(file, "%u %s",nHoles + MeshData::getInstance().getnbAddHoles(),"\n");
  unsigned iHole=0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned I=1; I<nShockPoints->at(ISH)-1; I++) {
    ++iHole;
    fprintf(file,"%u %s",iHole,"  ");
    fprintf(file,"%.16f %s %.16f %s",(*XYSh)(0,I,ISH),"  ",(*XYSh)(1,I,ISH),"\n");
   }
  }

  for (unsigned I=0; I<MeshData::getInstance().getnbAddHoles(); I++) {
   fprintf(file,"%u %s",iHole+I,"  ");
   for(unsigned j=0; j<2; j++) { fprintf(file,"%.16f %s",caddholes->at(j)," ");}
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::setAddress()
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

void WriteTriangle::setMeshData()
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

void WriteTriangle::setPhysicsData()
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
