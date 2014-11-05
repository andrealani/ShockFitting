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
  delete coor; delete zroe; delete celnod;
  delete celcel; delete edgptr; delete bndfac;
}

//--------------------------------------------------------------------------//

void ReadTriangle::setup()
{
  LogToScreen(VERBOSE, "ReadTriangle::setup() => start\n");

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

  setMeshData();
  setPhysicsData();

  if ((*firstRead)==1) {
   npoin->resize(2);
   nbfac->resize(2);
   nelem->resize(2);
   nedge->resize(2);
   nhole->resize(2);
  }

  logfile.Open(getClassName());

  ReadNode();

  ReadPoly();

  ReadEle();

  ReadNeigh();

  ReadEdge();

  logfile.Close();

  (*firstRead)=0;
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadNode()
{
  unsigned totsize,dim,iattr, states;
  unsigned idum;

  file.open(getNodeFile().c_str());
  file >> m_npoin >> dim >> states >> iattr;
  if (dim != (*ndim) || states > (*ndofmax)) {
   logfile("WARNING !!!!!");
   logfile("WARNING !!!!!");
   logfile("WARNING !!!!!");
   logfile("NDIM: ", dim, " NDOF: ", states, " in file ");
   logfile(getNodeFile().c_str(),"\n");
   logfile("NDIM: ", *ndof, " NDOF: ", *ndofmax, " in PhysicData\n");
   logfile("WARNING !!!!!");
   exit(1);
  }
   logfile("There are ", m_npoin, " gridpoints in ");
   logfile(getNodeFile().c_str(), "file\n");
   logfile("There are ", *ndof, " ndof in ");
   logfile(getNodeFile().c_str(), ".node file\n");
  if (iattr !=1) {
   cout << "ReadTriangle::error => <# of attributes> should be 1 in node file\n";
   exit(1);
  }

  if     ((*firstRead)==1)  { npoin->at(0) = m_npoin; }
  else if ((*firstRead)==0) { npoin->at(1) = m_npoin; }

  // resize vectors of the MeshData pattern
  // and assign dstarting pointers for arrays 2D
  // according to firstRead flag  
  totsize = 0;
  if      ((*firstRead)==1) {
   
   totsize0 = npoin->at(0) + 2 * (*nshmax) * (*npshmax);
   nodcod->resize(totsize0);
   coorVect->resize((*ndim) * totsize0);
   zroeVect->resize((*ndofmax) * totsize0);
   coor = new Array2D<double>( (*ndim) , totsize0, &coorVect->at(0) );
   zroe = new Array2D<double>( (*ndofmax) , totsize0, &zroeVect->at(0) );
  }
  else if  ((*firstRead)==0) {
   
  totsize = npoin->at(0) + npoin->at(1) + 4 * (*nshmax) * (*npshmax);
  zroeVect->resize((*ndofmax) * totsize);
  coorVect->resize((*ndim) * totsize);
  start = (*ndim) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  coor = new Array2D <double> ((*ndim),
                               (npoin->at(1) + 2 * (*nshmax) * (*npshmax)),
                               &coorVect->at(start));
  start = (*ndofmax) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  zroe = new Array2D <double>((*ndofmax),
                              (npoin->at(1) + 2 * (*nshmax) * (*npshmax)),
                              &zroeVect->at(start));
  }

  // fill zroe, coor and nodcod
  for (unsigned IPOIN=0; IPOIN < m_npoin; IPOIN++) {
   file >> idum;
  for(unsigned IA=0; IA < (*ndim); IA++) {
   file >> (*coor)(IA,IPOIN);
  }
  for(unsigned IA=0; IA < (*ndof); IA++) {
   file >> (*zroe)(IA,IPOIN);
  }
  if      ((*firstRead)==1)       { file >> nodcod->at(IPOIN); }
  else if ((*firstRead)==0)       { start = npoin->at(0) + 2 * (*nshmax) * (*npshmax);
                                    file >> nodcod->at(IPOIN+start); }
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
   cout << "ReadTriangle::error => node list should be empty in the ";
   cout << getPolyFile().c_str() << " file" << endl;
   exit(1);
  }
  file >> m_nbfac >> iattr;
  if (iattr!=1) {
   cout << "ReadTriangle::error => <# of attributes> should be 1 in the ";
   cout << getPolyFile().c_str() << " file" << endl;
   exit(1);
  }

  if      ((*firstRead)==1)  { nbfac->at(0) = m_nbfac; }
  else if ((*firstRead)==0)  { nbfac->at(1) = m_nbfac; }

  // resize vectors of the MeshData pattern
  // and assign dstarting pointers for arrays 2D
  // according to firstRead flag  
  if      ((*firstRead)==1) {

   totsize0 = nbfac->at(0) + 2 * (*nshmax) * (*neshmax);
   bndfacVect->resize(3 * totsize0);
   bndfac = new Array2D<int> (3,totsize0,&bndfacVect->at(0));
  }

  else if ((*firstRead)==0) { 

   totsize = nbfac->at(0) + nbfac->at(1) + 4 * (*nshmax) * (*neshmax);
   bndfacVect->resize(3 * totsize);
   start = 3* (nbfac->at(0) + 2 * (*nshmax) * (*neshmax));
   bndfac = new Array2D<int> (3,(nbfac->at(1) + 2 * (*nshmax) * (*neshmax)),
                              &bndfacVect->at(start));
  }

  logfile("There are",m_nbfac," polylines in ",getPolyFile(),"\n");

  for (unsigned IFACE=0; IFACE < m_nbfac; IFACE++) {
   file >> idum;
   for (unsigned IA=0; IA < 3; IA++) {file >> (*bndfac)(IA,IFACE);}
  }
  file >> m_nhole;
  if      ((*firstRead)==1) { nhole->at(0) = m_nhole; }
  else if ((*firstRead)==0) { nhole->at(1) = m_nhole; }

  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadEle()
{
  unsigned states, idum;

  file.open(getEleFile().c_str());

  file >> m_nelem >> *nvt >> states;
  if((*nvt)!=3) {
   cout << "ReadTriangle::error => <nodes per triangle> MUST be 3 in " ;
   cout << getEleFile().c_str() << " file" << endl;
   exit(1);
  }
  if(states!=0) {
   cout << "ReadTriangle::error => <# of attributes> MUST be 0 in ";
   cout << getEleFile().c_str() << " file" << endl;
   exit(1);
  }

  if      ((*firstRead)==1) { nelem->at(0) = m_nelem; }
  else if ((*firstRead)==0) { nelem->at(1) = m_nelem; }

  // resize vectors of the MeshData pattern
  // and assign dstarting pointers for arrays 2D
  // according to firstRead flag  
  if      ((*firstRead)==1) {

   totsize0 = nelem->at(0);
   celnodVect->resize((*nvt) * totsize0);
   celnod = new Array2D<int> ((*nvt), totsize0, &celnodVect->at(0));
  }

  else if ((*firstRead)==0) {

  totsize = nelem->at(0) + nelem->at(1);
  celnodVect->resize((*nvt) * totsize);
  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  }

  logfile("There are ",m_nelem," triangles in ",getEleFile(),"\n");

  //read data from .ele file and fill mesh array
  for (unsigned IELEM=0; IELEM < m_nelem; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA < *nvt; IA++) {file >> (*celnod)(IA,IELEM);}
  }
  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadNeigh()
{
  unsigned idum;

  file.open(getNeighFile().c_str());

  file >> m_nelem >> *nvt;
  if((*nvt) !=3) {
   cout << "ReadTriangle::error => <nodes per triangle> MUST be 3 in .neigh file \n";
   exit(1);
  }

  // resize vectors of the MeshData pattern
  // and assign dstarting pointers for arrays 2D
  // according to firstRead flag  
  if      ((*firstRead)==1) { 

   totsize0 = nelem->at(0);
   celcelVect->resize((*nvt)*totsize0);
   celcel = new Array2D<int> ((*nvt), totsize0, &celcelVect->at(0));
  }

  else if ((*firstRead)==0) {

   totsize = nelem->at(0) + nelem->at(1);
   celcelVect->resize((*nvt) * totsize);
   start = (*nvt) * nelem->at(0);
   celcel = new Array2D<int> ((*nvt), nelem->at(1), &celcelVect->at(start));
  }

  logfile("There are ",m_nelem," triangles in ",getNeighFile(),"\n");

  //read data from .neigh file and fill mesh array
  for (unsigned IELEM=0; IELEM<m_nelem; IELEM++) {
   file >> idum;
   for (unsigned IA=0; IA<(*nvt); IA++) { file >> (*celcel)(IA,IELEM); }
 }

  file.close();
}

//--------------------------------------------------------------------------//

void ReadTriangle::ReadEdge()
{
  unsigned iattr, idum;

  file.open(getEdgeFile().c_str());

  file >> m_nedge >> iattr;

  logfile("There are ",m_nedge," edge in ",getEdgeFile(),"\n");

  if (iattr!=1) {
   cout << "ReadTriangle::error => <# of attributes should be 1 in edge file>\n";
   exit(1);
  }

  if      ((*firstRead)==1)  { nedge->at(0) = m_nedge; }
  else if ((*firstRead)==0)  { nedge->at(1) = m_nedge; }

  // resize vectors of the MeshData pattern
  // and assign dstarting pointers for arrays 2D
  // according to firstRead flag  
  if      ((*firstRead)==1) {

   totsize0 = nedge->at(0);
   edgptrVect->resize(3 * totsize0);
   edgptr = new Array2D<int> (3, totsize0, &edgptrVect->at(0));
  }

  else if ((*firstRead)==0) {

   totsize = nedge->at(0) + nedge->at(1);
   edgptrVect->resize(3 * totsize);
   start = nedge->at(0) * 3;
   edgptr = new Array2D<int> (3, nedge->at(1), &edgptrVect->at(start));
  }

  //read data from .edge file and fill mesh array
  for (unsigned IFACE=0; IFACE<(m_nedge); IFACE++) {
   file >> idum;
    for (unsigned IA = 0; IA<3 ; IA++) {file >> (*edgptr)(IA,IFACE);}
  }
  file.close();
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getInputFiles() const
{
  if((*firstRead)==1) { assert(m_inputFile.size()==1);
                        fname->at(0) = m_inputFile[0]; }
  return fname->at(0);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getNodeFile()const
{
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="node") {
    if ((*firstRead)==0) { return string(getInputFiles()+".1."+m_fileTypes[i]);}
    else                 { return string(getInputFiles()+"."+m_fileTypes[i]);}
   }
  }
  LogToScreen(INFO, "ReadTriangle::error => in input.case no <node> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getPolyFile() const
{
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="poly") {
    if ((*firstRead)==0) { return string(getInputFiles()+".1."+m_fileTypes[i]);}
    else                 { return string(getInputFiles()+"."+m_fileTypes[i]);}
   }
  }
  LogToScreen(INFO, "ReadTriangle::error => in input.case no <poly> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getEleFile() const
{
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="ele") {
    if ((*firstRead)==0) { return string(getInputFiles()+".1."+m_fileTypes[i]);}
    else                 { return string(getInputFiles()+"."+m_fileTypes[i]);}
   }
  }
  LogToScreen(INFO, "ReadTriangle::error => in input.case no <ele> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getNeighFile() const
{
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
   if (m_fileTypes[i]=="neigh") {
    if ((*firstRead)==0) { return string(getInputFiles()+".1."+m_fileTypes[i]);}
    else                 { return string(getInputFiles()+"."+m_fileTypes[i]);}
   }
  }
  LogToScreen(INFO, "ReadTriangle::error => in input.case no <neigh> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

std::string ReadTriangle::getEdgeFile() const
{
  assert(m_fileTypes.size()==5);
  for (unsigned i=0; i<m_fileTypes.size(); i++) {
  if (m_fileTypes[i]=="edge") {
    if ((*firstRead)==0) { return string(getInputFiles()+".1."+m_fileTypes[i]);}
    else                 { return string(getInputFiles()+"."+m_fileTypes[i]);}
   }
  }
  LogToScreen(INFO,"ReadTriangle::error => in input.case no <edge> type found");
  exit(1);
}

//--------------------------------------------------------------------------//

void ReadTriangle::setMeshData ()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nedge = MeshData::getInstance().getData <vector<unsigned> > ("NEDGE");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nhole = MeshData::getInstance().getData <vector<unsigned> > ("NHOLE");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  celcelVect = MeshData::getInstance().getData <vector<int> >("CELCEL");
  edgptrVect = MeshData::getInstance().getData <vector<int> >("EDGPTR");
  fname = MeshData::getInstance().getData <vector<string> >("FNAME");
  firstRead = MeshData::getInstance().getData <unsigned> ("FirstRead");
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
