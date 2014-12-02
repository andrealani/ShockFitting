// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/CFmesh2Triangle.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Jcycl.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/Factory.hh"
#include "SConfig/ConfigFileReader.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CFmesh2Triangle, Converter>
CFmesh2TriangleProv("CFmesh2Triangle");

//----------------------------------------------------------------------------//

CFmesh2Triangle::CFmesh2Triangle(const std::string& objectName) :
  Converter(objectName)
{
  m_prim2param.name() = "DummyVariableTransformer";
}

//----------------------------------------------------------------------------//

CFmesh2Triangle::~CFmesh2Triangle()
{
  delete XY; delete zroe; delete celnod;
  delete celcel; delete bndfac;
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setup()
{
  LogToScreen(VERBOSE, "CFmesh2Triangle::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "CFmesh2Triangle::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::unsetup()
{
  LogToScreen(VERBOSE, "CFmesh2Triangle::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::configure(OptionMap& cmap, const std::string& prefix)
{
  Converter::configure(cmap, prefix);

  // assign strings on input.case file to variable transformer object
  m_prim2param.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;


  if (ConfigFileReader::isFirstConfig()) {
   m_prim2param.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                             getProvider(m_prim2param.name())
                             ->create(m_prim2param.name()));
  }

  // configure variable transformer object
  configureDeps(cmap, m_prim2param.ptr().get());
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::convert()
{
  LogToScreen (INFO, "CFmesh2Triangle::convert()\n");

  setMeshData();
  setPhysicsData();

  // read CFmesh format file
  LogToScreen(DEBUG_MIN, "CFmesh2Triangle::reading CFmesh format\n");
  readCFmeshFmt();

  // make the transformation fro primitive variables to dimensional 
  // Roe parameter vector variables
  m_prim2param.ptr()->transform(); 

  if (MeshData::getInstance().getVersion()=="original") {
   // write Triangle format files
   LogToScreen(DEBUG_MIN, "CFmesh2Triangle::writing Triangle format\n");
   writeTriangleFmt(); 
  }
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::readCFmeshFmt()
{
  // dummy variables
  unsigned LSKIP=0; unsigned ISKIP=0; unsigned value; 
  string skipver = "dummy"; string dummy; string strvalue;

  // list states values
  unsigned LIST_STATES;

  // boundary colours
  unsigned NCLR;

  // space dimensions read from CFmesh file
  unsigned ndim;

  // vector of GEOM_ENTS
  vector<unsigned> nFacB;

  // reading file
  ifstream file;


  string cfoutCFmesh = "cfout.CFmesh";

  file.open(string(cfoutCFmesh).c_str());
  file.precision(16);

  // read number of the first dummy strings
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open(string(cfoutCFmesh).c_str());
  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  while(ISKIP<(LSKIP)) { file >> skipver >> dummy;
                          ++ISKIP;      }
  file >> dummy >> ndim;                              // read !NB_DIM    ndim
  file >> dummy >> (*ndof);                           // read !NB_EQ     ndof
  file >> dummy >> npoin->at(1) >> dummy;             // read !NB_NODES  npoin 0
  file >> dummy >> dummy >> dummy;                    // read !NB_STATES nbstates 0
  file >> dummy >> nelem->at(1);                      // read !NB_ELEM   nelem
  file >> dummy >> dummy;                             // read !NB_ELEM_TYPES nbelemTypes

  PhysicsInfo::setnbDim(ndim);

  // read !GEOM_POLYORDER   geomPolyOrder
  // read !SOL_POLYORDER    solPolyOrder
  // read !ELEM_TYPES       Triag
  // read !NB_ELEM_PER_TYPE nbElemPerType
  // read !NB_NODES_PER_TYPE nbNodesPerType
  // read !NB_STATES_PER_TYPE nbStatesPerType
  unsigned I=0;
  while(I<6) { file >> dummy >> dummy;
               ++I;                    }
  // read !LIST_ELEM
  file >> dummy;
  // read elements of LIST_ELEM 
  I=0;
  while(I<nelem->at(1)) { 
   file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
   ++I;
  }

  // read !NB_TRSs nbTRSs
  file >> dummy >> NCLR;

  namebnd.resize(NCLR);
  nFacB.resize(NCLR);
  for(unsigned IFACE=0; IFACE<NCLR; IFACE++) {
   file >> dummy >> strvalue;   // read !TRS_NAME trsName
   namebnd.at(IFACE)= strvalue;
   file >> dummy >> dummy;      // read !NB_TRs nbTRs   
   file >> dummy >> value;      // read !NB_GEOM_ENTS nbgeomEnts
   nFacB.at(IFACE) = value;
   file >> dummy >> dummy;      // read !GEOM_TYPE geomType
   file >> dummy;               // read !LIST_GEOM_ENT
   for(unsigned IBND=0; IBND<nFacB.at(IFACE); IBND++) { 
    for(unsigned i=0; i<6; i++) {file >> dummy; }
   }
  }

  file.close();
    
  nbfac->at(1) = 0;
  for(unsigned IFACE=0; IFACE<NCLR; IFACE++) {
   nbfac->at(1) = nFacB.at(IFACE) + nbfac->at(1);
  }

  (*nvt) = PhysicsInfo::getnbDim() + 1;
  nhole->at(1) = 0; 

  // resize vectors of MeshData pattern with the new values read on CFmesh file
  // and assign starting pointers for arrays 2D
  resizeVectors();

  // read mesh points coordinates, mesh points status, the mesh connectivity,
  // the boundary structure and the solution from CFmesh file
  file.open(string(cfoutCFmesh).c_str());

  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  LSKIP = 0;
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");


  file.open(string(cfoutCFmesh).c_str());
  ISKIP=0;
  // read !NB_DIM ndim
  // read !NB_EQ ndof
  while(ISKIP<(2+LSKIP)) { file >> dummy >> dummy;
                            ++ISKIP;            }
  // read !NB_NODES
  // read !NB_STATES
  while(ISKIP<(4+LSKIP)) { file >> dummy >> dummy >> dummy;
                            ++ISKIP;  }
  // read !NB_ELEM
  // read !NB_ELEM_TYPES
  // read !GEOM_POLYORDER
  // read !SOL_POLYORDER
  // read !ELEM_TYPES
  // read !NB_ELEM_PER_TYPE 
  // read !NB_NODES_PER_TYPE
  // read !NB_STATES_PER_TYPE
  while(ISKIP<(12+LSKIP)) { file >> dummy >> dummy;
                            ++ISKIP;            }

  // read !LIST_ELEM
  file >> dummy;

  // read celnod values on !LIST_ELEM
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   for(unsigned IVERT=0; IVERT<(*nvt); IVERT++) {
    file >> value;
    (*celnod)(IVERT,IELEM) = value+1;
   }
   for(unsigned IVERT=0; IVERT<(*nvt); IVERT++) {
    file >> value; 
   }
  }

  file >> dummy >> dummy; // !NB_TRs numberTRs

  unsigned IFACE=0;
  double bndfac0, bndfac1;
  for(unsigned IC=0; IC<NCLR; IC++) { 
   ISKIP=1; 
   // read !TRS_NAME trsname
   // read !NB_TRs nbTrs
   // read !NB_GEOM_ENTS nbGeomEnts
   // read !GEOM_TYPE geomType
   while(ISKIP<5) { file >> dummy >> dummy;
                    ++ISKIP; } 
   // read !LIST_GEOM_ENT
   file >> dummy;

   // read geometry entities list that is boundary data
   for(unsigned IB=0; IB<nFacB.at(IC);IB++) {
    file >> dummy >> dummy >> bndfac0 >> bndfac1 >> dummy >> dummy;
    (*bndfac)(0,IFACE) = bndfac0 + 1;
    (*bndfac)(1,IFACE) = bndfac1 + 1;
    (*bndfac)(2,IFACE) = IC+1; // c++ indeces start from 0
    ++IFACE;
   }
  }

  // read !LIST_NODE and other dummy strings before the nodal coordinates
  file >> dummy; // read !EXTRA_VARS
  file >> dummy; // read !LIST_NODES
  if(dummy!="!LIST_NODE") { file >> dummy; }

  // read the nodal coordinates
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   file >> (*XY)(0,IPOIN) >> (*XY)(1,IPOIN); }

  file >> dummy >> LIST_STATES;
  if      (LIST_STATES==0) {
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned I=0; I<(*ndof); I++) { (*zroe)(I,IPOIN) = 0; }
   }
  }
  else if (LIST_STATES==1) {
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned I=0; I<(*ndof); I++) { file >> (*zroe)(I,IPOIN); }
   }
  }
  else { cout << "CFmesh2Triangle::error => !LIST_STATES must be 0 or 1 ";
         cout << "in cfout.CFmesh file\n"; }

  // fill nodcod vector
  setNodcod();

  // count number of boundary points
  countnbBoundaryNodes();

  file.close();
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::resizeVectors()
{
  totsize = npoin->at(0) + npoin->at(1) + 
            4 * PhysicsInfo::getnbShMax() *
                PhysicsInfo::getnbShPointsMax();
  nodcod->resize(totsize);
  unsigned startNodcod = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                            PhysicsInfo::getnbShPointsMax();
  // initialize nodcod of the shocked mesh
  for(unsigned I=startNodcod; I<totsize; I++) {
   nodcod->at(I) = 0; }

  zroeVect->resize(PhysicsInfo::getnbDofMax() * totsize);
  coorVect->resize(PhysicsInfo::getnbDim() * totsize);

  start = PhysicsInfo::getnbDim() * 
          (npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax());
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 * 
                             PhysicsInfo::getnbShMax() * 
                             PhysicsInfo::getnbShPointsMax()),
                             &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() *
          (npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double>(PhysicsInfo::getnbDofMax(),
                              (npoin->at(1)+2 * 
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                              &zroeVect->at(start));

  totsize = nelem->at(0) + nelem->at(1);
  celcelVect->resize((*nvt) * totsize);
  celnodVect->resize((*nvt) * totsize);

  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  celcel = new Array2D<int> ((*nvt), nelem->at(1), &celcelVect->at(start));

  totsize = nbfac->at(0) + nbfac->at(1) + 
            4 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax();
  bndfacVect->resize(3 * totsize);

  start = 3 * (nbfac->at(0) + 
            2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax());
  bndfac = new Array2D<int> (3,(nbfac->at(1) + 2 * PhysicsInfo::getnbShMax() *
                                PhysicsInfo::getnbShEdgesMax()),
                             &bndfacVect->at(start));
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setNodcod()
{
  int IPOIN, help;

  // write new nodcode values on the vector of the shocked mesh
  // in the fortran version this new vector is referred to index 1 (NODCOD(1))
  // here it is pushed back to the nodcod of the background mesh
  unsigned startNodcod = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                            PhysicsInfo::getnbShPointsMax();

  for(unsigned IFACE=0; IFACE<nbfac->at(1); IFACE++) {
   for(unsigned I=0; I<(*nvt)-1; I++) {
    IPOIN = (*bndfac)(I,IFACE);
    help = nodcod->at(startNodcod+IPOIN-1);
    nodcod->at(startNodcod+IPOIN-1) = help+1;
   }
  } 
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::countnbBoundaryNodes()
{
  // write new nodcode values on the vector of the shocked mesh
  // in the fortran version this new vector is referred to index 1 (NODCOD(1))
  // here it is pushed back to the nodcod of the background mesh
  unsigned startNodcod = 
   npoin->at(0) + 2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();

  unsigned nbBoundaryNodes=0;
  for(unsigned IPOIN=0;IPOIN<npoin->at(1); IPOIN++) {
   if(nodcod->at(startNodcod+IPOIN)!=0) { ++nbBoundaryNodes;}
  }
  nbpoin->at(1) = nbBoundaryNodes;
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::writeTriangleFmt()
{
  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  string NBND;
 
  // writing file (Triangle node file to be overwritten)
  FILE* file;

  // take new nodcode values from the new nodcod vector of the shocked mesh
  // in the fortran version this new vector is referred to index 1 (NODCOD(1))
  // here it is pushed back to the nodcod of the background mesh
  unsigned startNodcod = 
   npoin->at(0) +2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();

  // write on .node file
  dummystring = fname->str()+".1.node";
  file = fopen(dummystring.c_str(),"w");

  fprintf(file,"%u %s %u",npoin->at(1)," ",PhysicsInfo::getnbDim());
  fprintf(file,"%s %u %s"," ",(*ndof)," 1\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   fprintf(file,"%u %s",IPOIN+1," ");
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++)
    { fprintf(file,"%.16f %s",(*XY)(IA,IPOIN)," "); }
   for(unsigned IA=0; IA<(*ndof); IA++) 
    { fprintf(file,"%.16f %s",(*zroe)(IA,IPOIN)," "); }
   fprintf(file,"%u %s",nodcod->at(startNodcod+IPOIN),"\n");
  }

  fclose(file);

  // write on .poly file
  dummystring = fname->str()+".1.poly";
  file = fopen(dummystring.c_str(),"w");

  fprintf(file,"%s %u %s", "0 ", PhysicsInfo::getnbDim(), " 0 1\n");
  fprintf(file,"%u %s",nbfac->at(1)," 1\n");
  for(unsigned IFACE=0; IFACE<nbfac->at(1); IFACE++) {
   NBND = namebnd.at((*bndfac)(2,IFACE)-1); // c++ indeces start from 0
   if(NBND=="InnerSup" || NBND=="InnerSub") { NBND="10"; }
   fprintf(file,"%u %s",IFACE+1," ");
   fprintf(file,"%i %s %i",(*bndfac)(0,IFACE)," ",(*bndfac)(1,IFACE));
   fprintf(file,"%s %s %s"," ",NBND.c_str()," \n");
  }

  fprintf(file,"%s","0\n"); // write number of holes   
  fclose(file);
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbpoin = MeshData::getInstance().getData <vector<unsigned> > ("NBPOIN");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nhole = MeshData::getInstance().getData <vector<unsigned> > ("NHOLE");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  celcelVect = MeshData::getInstance().getData <vector<int> >("CELCEL");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  fname = MeshData::getInstance().getData <stringstream>("FNAME");
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

