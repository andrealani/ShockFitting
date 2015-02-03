// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteBackTriangle.hh"
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
ObjectProvider<WriteBackTriangle, WritingMesh> 
  writeBackTriangleProv("WriteBackTriangle");

//--------------------------------------------------------------------------//

WriteBackTriangle::WriteBackTriangle(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteBackTriangle::~WriteBackTriangle()
{
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setup()
{
  LogToScreen(VERBOSE,"WriteBackTriangle::setup() => start\n");

  LogToScreen(VERBOSE,"WriteBackTriangle::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::unsetup()
{
  LogToScreen(VERBOSE,"WriteBackTriangle::unsetup()\n");
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::write()
{
  LogToScreen(INFO,"WriteBackTriangle::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  *fnameBack = "na99";
  string dummyfile = *fnameBack+".node";

  file = fopen(dummyfile.c_str(),"w");

  unsigned ilist = npoin->at(0);
  unsigned iPoin;

  fprintf(file, "%u %s %u %s %u %s",ilist,"  ",PhysicsInfo::getnbDim(),"  ",*ndof,"  1\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   iPoin = IPOIN+1;
   if      (nodcod->at(IPOIN)==-1) {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XY)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) { 
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%s","0\n");
   }
   else if (nodcod->at(IPOIN)==-2) {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XY)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {  
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%s","2\n");
   }
   else {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XY)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {  
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%u %s",nodcod->at(IPOIN), "\n");
   }
  }
  
  fclose(file);  

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setAddress()
{
  Zroe = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                             npoin->at(0),
                             &zroeVect->at(0));
  XY = new Array2D<double>(PhysicsInfo::getnbDim(), npoin->at(0),
                           &coorVect->at(0));
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::freeArray()
{
  delete Zroe; delete XY;
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  fnameBack = MeshData::getInstance().getData <string>("FNAMEBACK");
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
