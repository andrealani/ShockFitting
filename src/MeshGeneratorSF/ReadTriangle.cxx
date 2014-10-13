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
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/FileLogManip.hh"
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
  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "ReadTriangle::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReadTriangle::unsetup()
{
  LogToScreen(VERBOSE, "ReadTriangle::unsetup()\n");
  logfile.Close();
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
  unsigned totsize,dim,iattr, states;
  unsigned idum, h1, h2;

  file.open(getNodeFile().c_str());
  file >> *npoin >> dim >> states >> iattr;
  if (dim != (*ndim) || states > (*ndofmax)) {
   logfile("WARNING !!!!!");
   logfile("WARNING !!!!!");
   logfile("WARNING !!!!!");
   logfile("NDIM: ", dim, " NDOF: ", states, " in file *.node\n");
   logfile("NDIM: ", *ndof, " NDOF: ", *ndofmax, " in PhysicData.hh\n");
   logfile("WARNING !!!!!");
   exit(1);
  }
   logfile("There are ", *npoin, " gridpoints in ", ".node file\n"); // fwork!!
   logfile("There are ", *ndof, " ndof in ", ".node file\n"); // fwork!!
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

  logfile("There are",*nbfac," polylines in ",getPolyFile(),"\n");

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
  unsigned nvt, states, idum;

  file.open(getEleFile().c_str());

  file >> *nelem >> nvt >> states;
  if(nvt!=3) {
   cout << "<nodes per triangle> MUST be 3 in " << getEleFile().c_str() << " file" << endl;
   exit(1);
  }
  if(states!=0) {
   cout << "<# of attributes> MUST be 0 in " << getEleFile().c_str() << " file" << endl;
   exit(1);
  }

  // resize array with correct size value
  celnod->resize(nvt,*nelem);

  logfile("There are ",*nelem," triangles in ",getEleFile(),"\n");

  //read data from .ele file and fill mesh array
  for (unsigned IELEM=0; IELEM < *nelem; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA < nvt; IA++) {file >> (*celnod)(IA,IELEM);}
  }
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadNeigh()
{
  unsigned nvt, idum;

  file.open(getNeighFile().c_str());

  file >> *nelem >> nvt;
  if(nvt !=3) {
   cout << "<nodes per triangle> MUST be 3 in .neigh file " << endl;
   exit(1);
  }

  //resize array with correct size value
  celcel->resize(nvt,*nelem);

  logfile("There are ",*nelem," triangles in ",getNeighFile(),"\n");

  //read data from .neigh file and fill mesh array
 for (unsigned IELEM=0; IELEM<*nelem; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA<nvt; IA++) {file >> (*celcel)(IA,IELEM);}
 }
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadEdge()
{
  unsigned iattr, idum;

  file.open(getEdgeFile().c_str());

  file >> *nedge >> iattr;

  logfile("There are ",*nedge," edge in ",getEdgeFile(),"\n");

  if (iattr!=1) {
   cout << "<# of attributes should be 1 in edge file>" << endl;
   exit(1);
  }

  // resize array with correct size value
  edgptr->resize(3,(*nedge));

  //read data from .edge file and fill mesh array
  for (unsigned IFACE=0; IFACE<(*nedge); IFACE++) {
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
  nedge = MeshData::getInstance().getData <unsigned> ("NEDGE");
  nelem = MeshData::getInstance().getData <unsigned> ("NELEM");
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
