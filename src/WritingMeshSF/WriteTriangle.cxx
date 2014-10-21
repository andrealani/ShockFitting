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
  LogToScreen(VERBOSE, "WriteTriangle::unsetup() => start\n");

  LogToScreen(VERBOSE, "WriteTriangle::unsetup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteTriangle::write()
{
  LogToScreen(INFO, "WriteTriangle::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  // write node file
  file.open("na00.2.node");

  ilist = (*npoin) + 2 * (*nshmax) * (*npshmax);

  M02M1.resize(ilist);
  M12M0.resize(ilist);

  // set map vector for nodcod
  setMapVectorForNodcod();

  icount = *npoin;

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();


  ilist = TNPOIN;
  file << ilist << " " << *ndim << " " << *ndof << " 1\n";
  
  // write mesh points coordinates and status on Triangle file
  writeMeshVariables();

  icount = (*npoin);

  // write upstream shock points coordinates and status on Triangle file
  writeUpstreamStatus();

  // write downstream shock points coordinates and status on Triangle file
  writeDownstreamStatus();

  file.close();

  // write poly file
  file.open("na00.2.poly");

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
  for (unsigned IPOIN=0; IPOIN< (*npoin); IPOIN++) {
   if (nodcod->at(IPOIN)>=0) { M02M1.at(IPOIN) = TNPOIN;
                               M12M0.at(TNPOIN) = IPOIN; 
			       ++TNPOIN;                }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::setMapVectorForNodcodSh()
{
  for (unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    if ((*w_NodCodSh)(I,ISH)==10) { M02M1.at(icount) = TNPOIN;
                                    M12M0.at(TNPOIN) = icount;
                                    ++TNPOIN;                 }
    ++icount;
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeMeshVariables()
{
  unsigned h = 0;
  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   if(nodcod->at(IPOIN)>=0) {
    file << M02M1.at(IPOIN) << "\n";
    for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XY)(IA,IPOIN) << " ";}
    file << "\n";
    for(unsigned IA=0; IA<(*ndof); IA++) {
     file << zroe->at(h) << " "; h++; }
    file << "\n";
    file << nodcod->at(IPOIN);
   }
  file << "\n";
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeUpstreamStatus()
{
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    if((*w_NodCodSh)(I,ISH)==10) {
     file << M02M1.at(icount) << "\n";
     for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XYShu)(IA,I,ISH) << " ";}
     file << "\n";
     for(unsigned IA=0; IA<(*ndof); IA++) {
     file << (*w_ZRoeShu)(IA,I,ISH) << " ";}
     file << "\n";
     file << (*w_NodCodSh)(I,ISH);
    }
    ++icount;
    file << "\n";
   }
   file << "\n";
  }
  file << "\n";
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeDownstreamStatus()
{
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    if((*w_NodCodSh)(I,ISH)==10) {
     file << M02M1.at(icount) << "\n";
     for(unsigned IA=0; IA<(*ndim); IA++) {
     file << (*w_XYShd)(IA,I,ISH) << " ";}
     file << "\n";
     for(unsigned IA=0; IA<(*ndof); IA++) {
     file << (*w_ZRoeShd)(IA,I,ISH) << " ";}
     file << "\n";
     file << (*w_NodCodSh)(I,ISH);
    }
    icount++;
    file << "\n";
   }
   file << "\n";
  }
  file << "\n";
}

//--------------------------------------------------------------------------//

void WriteTriangle::writeBndfac()
{
  int IBC;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(2,IFACE);
   if(IBC>0) { ++ICHECK; }
  }
  file << ICHECK << " 1\n";

  ICHECK=0;

  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) { 
   IBC = (*bndfac)(2,IFACE);
   if(IBC>0) {
    ++ICHECK; 
    file << ICHECK  << " ";
    for(unsigned IA=0; IA<1; IA++) {
     file << M02M1.at((*bndfac)(IA,IFACE)) << " ";
    }
   file << IBC << "\n";
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::computenbHoles()
{
  for(unsigned ISH=0; ISH<(*w_nShocks); ISH++) {
   nHoles = nHoles + w_nShockPoints->at(ISH)-2;
  }

  file << nHoles + (*naddholes);

  unsigned iHole=0;
  for(unsigned ISH=0; ISH<(*w_nShocks); ISH++) {
   for(unsigned I=0; I<w_nShockPoints->at(ISH); I++) {
    ++iHole;
    file << iHole << " ";
    file << (*w_XYSh)(0,I,ISH) << " " << (*w_XYSh)(1,I,ISH) << "\n";
   }
  }

  for (unsigned I=0; I<(*naddholes); I++) {
   file << iHole+I << " ";
   for(unsigned j=0; j<2; j++) { file << caddholes->at(j) << " ";}
  }
}

//--------------------------------------------------------------------------//

void WriteTriangle::setAddress()
{
  unsigned start;
  start = 0;
  w_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
  start = (*npoin);
  w_NodCodSh = new Array2D <int> ((*npshmax),(*nshmax),&nodcod->at(start));
  start = (*npoin)*(*ndof);
  w_ZRoeShu = 
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
  start = (*npoin) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  w_ZRoeShd = 
    new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroe->at(start));
  start = (*npoin) * (*ndim);
  w_XYShu =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
  start = (*npoin) * (*ndim) + (*npshmax) * (*nshmax) * (*ndim);
  w_XYShd =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
}

//--------------------------------------------------------------------------//

void WriteTriangle::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  nbfacSh = MeshData::getInstance().getData<unsigned>("NBFACSH");
  naddholes = MeshData::getInstance().getData<unsigned>("Naddholes");
  caddholes = 
    MeshData::getInstance().getData <std::vector<double> > ("CADDholes");
  nodcod = MeshData::getInstance().getData< std::vector<int> >("NODCOD");
  zroe = MeshData::getInstance().getData< std::vector<double> >("ZROE");
  coor = MeshData::getInstance().getData< std::vector<double> >("COOR");
  bndfac = MeshData::getInstance().getData< Array2D<int> >("BNDFAC");
  celnod = MeshData::getInstance().getData< Array2D<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void WriteTriangle::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  w_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  w_nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  w_nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  w_XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
