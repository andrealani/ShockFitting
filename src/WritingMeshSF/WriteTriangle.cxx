// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteTriangle.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
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
  TNPOIN=0; ICHECK=0; nHoles=0;
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
  fname->at(0) = "na00001";

  // write node file
  dummyfile = fname->at(0)+".node";
  file.open(dummyfile.c_str());

  file.precision(20);

  ilist = npoin->at(0) + 2 * (*nshmax) * (*npshmax);

  M02M1->resize(ilist+1); // c++ indeces start from 0
  M12M0->resize(ilist+1); // c++ indeces start from 0

  // set map vector for nodcod
  setMapVectorForNodcod();

  icount = npoin->at(0);

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  ilist = TNPOIN;

  file << ilist << " " << *ndim << " " << *ndof << " 1\n";
  
  // write mesh points coordinates and status on Triangle file
  writeMeshVariables();

  icount = npoin->at(0);

  // write upstream shock points coordinates and status on Triangle file
  writeUpstreamStatus();

  // write downstream shock points coordinates and status on Triangle file
  writeDownstreamStatus();

  file.close();

  // write poly file
  dummyfile = fname->at(0)+".poly";
  file.open(dummyfile.c_str());

  file << "0 " << (*ndim) << " 0 1 \n";
  
  // set map vector for bndfac and write bndfac on poly file
  writeBndfac();

  // compute number of holes
  computenbHoles();

  file.close();
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
  for (unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    ++icount;

    if ((*w_NodCodSh)(I,ISH)==10) { ++TNPOIN;
                                    M02M1->at(icount) = TNPOIN;
                                    M12M0->at(TNPOIN) = icount;}
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeMeshVariables()
{
  unsigned h = 0;
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if(nodcod->at(IPOIN)>=0) {
    file << M02M1->at(IPOIN+1) << "  "; // c++ indeces start from 0
    for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XY)(IA,IPOIN) << "  ";}
    for(unsigned IA=0; IA<(*ndof); IA++) {
     file << zroeVect->at(h) << "  "; h++; }
    file << nodcod->at(IPOIN) << endl;
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeUpstreamStatus()
{
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    ++icount;
    if((*w_NodCodSh)(I,ISH)==10) {
     file << M02M1->at(icount) << "  ";
     for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XYShu)(IA,I,ISH) << "  ";}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     file << (*w_ZRoeShu)(IA,I,ISH) << "  ";}
     file << (*w_NodCodSh)(I,ISH) << endl;
    }
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeDownstreamStatus()
{
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
     icount++;
    if((*w_NodCodSh)(I,ISH)==10) {
     file << M02M1->at(icount) << "  ";
     for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XYShd)(IA,I,ISH) << "  ";}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     file << (*w_ZRoeShd)(IA,I,ISH) << "  ";}
     file << (*w_NodCodSh)(I,ISH) << endl;
    }
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeBndfac()
{
  int IBC;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(2,IFACE);
   if(IBC>0) { ++ICHECK; }
  }
  file << ICHECK << "  1" << endl;

  ICHECK=0;


  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) { 
   IBC = (*bndfac)(2,IFACE);



   if(IBC>0) {
    ++ICHECK; 
    file << ICHECK  << "  ";
    for(unsigned IA=0; IA<2; IA++) {
     file << M02M1->at((*bndfac)(IA,IFACE)) << "  ";
    }
   file << IBC << endl;
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::computenbHoles()
{
  for(unsigned ISH=0; ISH<(*w_nShocks); ISH++) {
   nHoles = nHoles + w_nShockPoints->at(ISH)-2;
  }

  file << nHoles + (*naddholes) << endl;

  unsigned iHole=0;
  for(unsigned ISH=0; ISH<(*w_nShocks); ISH++) {
   for(unsigned I=1; I<w_nShockPoints->at(ISH)-1; I++) {
    ++iHole;
    file << iHole << "  ";
    file << (*w_XYSh)(0,I,ISH) << "  " << (*w_XYSh)(1,I,ISH) << endl;
   }
  }

  for (unsigned I=0; I<(*naddholes); I++) {
   file << iHole+I << "  ";
   for(unsigned j=0; j<2; j++) { file << caddholes->at(j) << " ";}
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::setAddress()
{
  unsigned start; unsigned totsize;
  start = 0;
  w_XY = new Array2D <double> ((*ndim),npoin->at(0),&coorVect->at(start));
  totsize = nbfac->at(0) + 2 * (*nshmax) * (*neshmax);
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  celnod = new Array2D<int> ((*nvt), nelem->at(0), &celnodVect->at(start)); 
  start = npoin->at(0);
  w_NodCodSh = new Array2D <int> ((*npshmax),(*nshmax),&nodcod->at(start));
  start = npoin->at(0) * (*ndof);
  w_ZRoeShu = 
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroeVect->at(start));
  start = npoin->at(0) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  w_ZRoeShd = 
    new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroeVect->at(start));
  start = npoin->at(0) * (*ndim);
  w_XYShu =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coorVect->at(start));
  start = npoin->at(0) * (*ndim) + (*npshmax) * (*nshmax) * (*ndim);
  w_XYShd =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coorVect->at(start));
}

//--------------------------------------------------------------------------//

void WriteTriangle::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem =MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbfacSh = MeshData::getInstance().getData<unsigned>("NBFACSH");
  naddholes = MeshData::getInstance().getData<unsigned>("Naddholes");
  caddholes = 
    MeshData::getInstance().getData <vector<double> > ("CADDholes");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  fname = MeshData::getInstance().getData <vector <string> > ("FNAME");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
}

//--------------------------------------------------------------------------//

void WriteTriangle::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  neshmax = PhysicsData::getInstance().getData <unsigned> ("NESHMAX");
  w_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  w_nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  w_nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  w_XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
