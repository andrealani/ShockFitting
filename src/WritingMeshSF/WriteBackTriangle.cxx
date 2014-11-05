// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteBackTriangle.hh"
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

  fnameBack->at(0) = "na99";
  string dummyfile = fnameBack->at(0)+".node";
  file.open(dummyfile.c_str());
  file.precision(20);

  unsigned ilist = npoin->at(0);
  unsigned iPoin;

  file << ilist << "  " << *ndim << "  " << *ndof << "  1\n";
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   iPoin = IPOIN+1;
   if      (nodcod->at(IPOIN)==-1) {
    file << iPoin << " ";
    for(unsigned IA=0; IA<(*ndim); IA++) { file << (*XY)(IA,IPOIN) << " "; }
    for(unsigned IA=0; IA<(*ndof); IA++) { file << (*Zroe)(IA,IPOIN) << " "; }
    file << "0\n";
   }
   else if (nodcod->at(IPOIN)==-1) {
    file << iPoin << " ";
    for(unsigned IA=0; IA<(*ndim); IA++) { file << (*XY)(IA,IPOIN) << " "; }
    for(unsigned IA=0; IA<(*ndof); IA++) { file << (*Zroe)(IA,IPOIN) << " "; }
    file << "2\n";
   }
   else {
    file << iPoin << " ";
    for(unsigned IA=0; IA<(*ndim); IA++) { file << (*XY)(IA,IPOIN) << " "; }
    for(unsigned IA=0; IA<(*ndof); IA++) { file << (*Zroe)(IA,IPOIN) << " "; }
    file << nodcod->at(IPOIN) << "\n";
   }
  }
  
  file.close();  
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setAddress()
{
  Zroe = new Array2D<double>((*ndof),npoin->at(0),&zroeVect->at(0));
  XY = new Array2D<double>((*ndim),npoin->at(0),&coorVect->at(0));
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  fnameBack = MeshData::getInstance().getData <vector<string> >("FNAMEBACK");
}

//--------------------------------------------------------------------------//

void WriteBackTriangle::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
