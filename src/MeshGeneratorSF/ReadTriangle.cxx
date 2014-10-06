// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReadTriangle.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Common/MeshData.hh"
#include "Common/PhysicsData.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReadTriangle, MeshGenerator> readTriangleProv("ReadTriangle");

//--------------------------------------------------------------------------//

ReadTriangle::ReadTriangle(const std::string& objectName) :
  MeshGenerator(objectName)
{
  m_fileTypes = vector<string>();
  addOption("FileTypes",&m_fileTypes,
             "List of file types names");
}

//--------------------------------------------------------------------------//

ReadTriangle::~ReadTriangle()
{
}

//--------------------------------------------------------------------------//

void ReadTriangle::setup()
{
  LogToScreen(VERBOSE, "ReadTriangle::setup() => start\n");

  setMeshData();
  setPhysicsData();

  LogToScreen(VERBOSE, "ReadTriangle::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReadTriangle::unsetup()
{
  LogToScreen(VERBOSE, "ReadTriangle::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReadTriangle::generate()
{
  LogToScreen(INFO, "ReadTriangle::generate()\n");

  ReadNode();
  ReadPoly();
  ReadEle();
  ReadNeigh();
  ReadEdge();

}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadNode()
{
  unsigned totsize,dim,iattr;
  unsigned idum, h1, h2;

  file.open(getNodeFile().c_str());
  file >> *npoin >> dim >> *ndof >> iattr;
  if (dim != (*ndim) || (*ndof) > (*ndofmax)) {
//   readmesh << "WARNING !!!!!" << endl;
//   readmesh << "WARNING !!!!!" << endl;
//   readmesh << "NDIM: " << dim << " NDOF: " << states << " in file *.node" << endl;
//   readmesh << "NDIM: " << ndof << " NDOF: " << *ndofmax << " in PhysicData.hh" << endl;
//   readmesh << "WARNING !!!!!" << endl;
   exit(1);
  }
//  readmesh << "There are " << *npoin << " gridpoints int " << fwork << endl;
//  readmesh << "There are " << *ndof << " degrees of freedom in " << fwork << endl;
  if (iattr !=1) {
   cout << "# of attributes should be 1 in node file" << endl;
   exit(1);
  }
 
  //resize vectors with correct size values
  totsize = (*npoin) + 2 * (*nshmax) * (*npshmax);
  nodcod->resize(totsize);
  totsize = (*ndim) * ((*npoin) + 2 * (*nshmax) * (*npshmax));
  coor->resize(totsize);
  totsize = (*ndofmax) * ((*npoin) + 2 * (*nshmax) * (*npshmax));
  zroe->resize(totsize);

  // fill zroe, coor and nodcod
  h1=0; h2=0;
  for (unsigned IPOIN=0; IPOIN < (*npoin); IPOIN++) {
   file >> idum;
  for(unsigned IA=0; IA < (*ndim); IA++) {
   file >> coor->at(h1);h1++;
  }
  for(unsigned IA=0; IA < (*ndof); IA++) {
   file >> zroe->at(h2); h2++;
  }
  file >> nodcod->at(IPOIN);
 }
 file.close();

}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadPoly()
{
  unsigned totsize,ia,idum,idum1,idum2,iattr;

  file.open(getPolyFile().c_str());

  file >> ia >> idum >> idum1 >> idum2;
  if (ia!=0) {
   cout << "node list should be empty in the poly file" << endl;
   exit(1);
  }
  file >> *nbfac >> iattr;
  if (iattr!=1) {
   cout << "<# of attributes should be 1 in poly file>" << endl;
   exit(1);
  }

  // resize array with correct size values
  totsize = (*nbfac) + 2 * (*nshmax) * (*neshmax); // leave room for duplicated nodes
  bndfac->resize(3,totsize);

  //readmesh << "There are " << *nbfac << " polylines in " << getPolyFile().c_str() << "  file"  << endl;

  for (unsigned IFACE=0; IFACE < (*nbfac); IFACE++) {
   file >> idum;
   for (unsigned IA=0; IA < 3; IA++) {file >> (*bndfac)(IA,IFACE);}
  }
  file >> *nhole;
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadEle()
{
  unsigned NELEM, NVT, states, idum;

  file.open(getEleFile().c_str());

  file >> NELEM >> NVT >> states;
  if(NVT!=3) {
   cout << "<nodes per triangle> MUST be 3 in " << getEleFile().c_str() << " file" << endl;
   exit(1);
  }
  if(states!=0) {
   cout << "<# of attributes> MUST be 0 in " << getEleFile().c_str() << " file" << endl;
   exit(1);
  }

  // resize array with correct size value
  celnod->resize(NVT,NELEM);

//  readmesh << "There are " << NELEM << " triangles in " << getEleFile().c_str() << " file" << endl;

  //read data from .ele file and fill mesh array
  for (unsigned IELEM=0; IELEM < NELEM; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA < NVT; IA++) {file >> (*celnod)(IA,IELEM);}
  }
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadNeigh()
{
  unsigned NELEM, NVT, idum;

  file.open(getNeighFile().c_str());

  file >> NELEM >> NVT;
  if(NVT !=3) {
   cout << "<nodes per triangle> MUST be 3 in .neigh file " << endl;
   exit(1);
  }

  //resize array with correct size value
  celcel->resize(NVT,NELEM);

// readmesh << "There are " << NELEM << " triangles in " << getNeighFile().c_str() << " file" << endl;

  //read data from .neigh file and fill mesh array
 for (unsigned IELEM=0; IELEM<NELEM; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA<NVT; IA++) {file >> (*celcel)(IA,IELEM);}
 }
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadEdge()
{
  unsigned NEDGE, iattr, idum;

  file.open(getEdgeFile().c_str());

  file >> NEDGE >> iattr;
 // readmesh << "There are " << NEDGE << " edge in " << getEdgeFile().c_str() << " file" << endl;
  if (iattr!=1) {
   cout << "<# of attributes should be 1 in edge file>" << endl;
   exit(1);
  }

  // resize array with correct size value
  edgptr->resize(3,NEDGE);

  //read data from .edge file and fill mesh array
  for (unsigned IFACE=0; IFACE<NEDGE; IFACE++) {
   file >> idum;
    for (unsigned IA = 0; IA<3 ; IA++) {file >> (*edgptr)(IA,IFACE);}
  }
  file.close();
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getNodeFile()const
{
  using namespace std;
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="node") {return string(getInputFiles()+"."+m_fileTypes[i]);}
  }
  LogToScreen(INFO, "ReadTriangle::ReadNode => ERROR: no <.node> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getPolyFile() const
{
  using namespace std;
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="poly") {return string(getInputFiles()+"."+m_fileTypes[i]);}
  }
  LogToScreen(INFO, "ReadTriangle::ReadPoly => ERROR: no <.poly> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getEleFile() const
{
  using namespace std;
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="ele") {return string(getInputFiles()+"."+m_fileTypes[i]);}
  }
  LogToScreen(INFO, "ReadTriangle::ReadEle => ERROR: no <.ele> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getNeighFile() const
{
  using namespace std;
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="neigh") {return string(getInputFiles()+"."+m_fileTypes[i]);}
  }
  LogToScreen(INFO, "ReadTriangle::ReadNeigh => ERROR: no <.neigh> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getEdgeFile() const
{
  using namespace std;
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
  if (m_fileTypes[i]=="edge") {return string(getInputFiles()+"."+m_fileTypes[i]);}
  }
  return string("ReadTriangle::ReadEdge => ERROR: no <.edge> type found");
}

//--------------------------------------------------------------------------//

void ReadTriangle::setMeshData ()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  nbpoin = MeshData::getInstance().getData <unsigned> ("NBPOIN");
  nhole = MeshData::getInstance().getData <unsigned> ("NHOLE");
  nodcod = MeshData::getInstance().getData< std::vector<int> >("NODCOD");
  zroe = MeshData::getInstance().getData< std::vector<double> >("ZROE");
  coor = MeshData::getInstance().getData< std::vector<double> >("COOR");
  bndfac = MeshData::getInstance().getData< Array2D<int> >("BNDFAC");
  celnod = MeshData::getInstance().getData< Array2D<int> >("CELNOD");
  celcel = MeshData::getInstance().getData< Array2D<int> >("CELCEL");
  edgptr = MeshData::getInstance().getData< Array2D<int> >("EDGPTR");
}

//--------------------------------------------------------------------------//

void ReadTriangle::setPhysicsData ()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  neshmax = PhysicsData::getInstance().getData <unsigned> ("NESHMAX");
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
