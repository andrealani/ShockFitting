// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/CFmesh2Triangle.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"
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

  setAddress();

  // read CFmesh format file
  LogToScreen(DEBUG_MIN, "CFmesh2Triangle::reading CFmesh format\n");
  readCFmeshFmt();

  // make the transformation fro primitive variables to dimensional 
  // Roe parameter vector variables
  m_prim2param.ptr()->transform(); 

  // write Triangle format files
  LogToScreen(DEBUG_MIN, "CFmesh2Triangle::writing Triangle format\n");
  writeTriangleFmt(); 
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

  // vector of GEOM_ENTS
  vector<unsigned> nFacB;

  // reading file
  ifstream file;
  file.open("cfout.CFmesh_prova");

  // read number of the first dummy strings
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open("cfout.CFmesh_prova");
  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  while(ISKIP<(LSKIP)) { file >> skipver >> dummy;
                          ++ISKIP;      }
  file >> dummy >> (*ndim);           // read !NB_DIM    ndim
  file >> dummy >> (*ndof);           // read !NB_EQ     ndof
  file >> dummy >> (*npoin) >> dummy; // read !NB_NODES  npoin 0
  file >> dummy >> dummy >> dummy;    // read !NB_STATES nbstates 0
  file >> dummy >> (*nelem);          // read !NB_ELEM   nelem
  file >> dummy >> dummy;             // read !NB_ELEM_TYPES nbelemTypes
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
  while(I<(*nelem)) { 
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
    
  (*nbfac) = 0;
  for(unsigned IFACE=0; IFACE<NCLR; IFACE++) {
   (*nbfac) = nFacB.at(IFACE) + (*nbfac);
  }

  (*nvt) = (*ndim) + 1;
  (*nHole) = 0; 

  // resize vector and array with the new values read on CFmesh file
  resizeArrays();

  // set address with the new values read on CFmesh files for mesh points
  setAddress();

  // read mesh points coordinates, mesh points status, the mesh connectivity,
  // the boundary structure and the solution from CFmesh file
  file.open("cfout.CFmesh_prova");

  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  LSKIP = 0;
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open("cfout.CFmesh_prova");
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
  for(unsigned IELEM=0; IELEM<(*nelem); IELEM++) {
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
  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   file >> (*c_XY)(0,IPOIN) >> (*c_XY)(1,IPOIN); }

  file >> dummy >> LIST_STATES;
  if      (LIST_STATES==0) {
   for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
    for(unsigned I=0; I<(*ndof); I++) { (*c_Zroe)(I,IPOIN) = 0; }
   }
  }
  else if (LIST_STATES==1) {
   for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
    for(unsigned I=0; I<(*ndof); I++) { file >> (*c_Zroe)(I,IPOIN); }
   }
  }
  else { cout << "CFmesh2Triangle::error => !LIST_STATES must be 0 or 1 ";
         cout << "in cfout.CFmesh file\n"; }

  // fill nodcod vector
  setNodcod();

  // count number of boundary points
  countnbBoundaryNodes();

}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::resizeArrays()
{
  unsigned totsize;
  totsize = (*npoin) * (*ndim);
  coor->resize(totsize);
  totsize = (*npoin) * (*ndof);
  zroe->resize(totsize);
  totsize = (*npoin);
  nodcod->resize(totsize);
  celnod->resize((*nvt), (*nelem));
  celcel->resize((*nvt), (*nelem));
  bndfac->resize(3, (*nbfac));
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setNodcod()
{
  int IPOIN, help;

  for(unsigned i=0; i<(*npoin); i++) { nodcod->at(i)=0;}

  for(unsigned IFACE=0; IFACE<(*nbfac); IFACE++) {
   for(unsigned I=0; I<(*nvt)-1; I++) {
    IPOIN = (*bndfac)(I,IFACE);
    help = nodcod->at(IPOIN-1);
    nodcod->at(IPOIN-1) = help+1;
   }
  } 
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::countnbBoundaryNodes()
{
  unsigned nbBoundaryNodes=0;
  for(unsigned IPOIN=0;IPOIN<(*npoin); IPOIN++) {
   if(nodcod->at(IPOIN)!=0) { ++nbBoundaryNodes;}
  }
  (*nbpoin) = nbBoundaryNodes;
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::writeTriangleFmt()
{
  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  string NBND;
 
  // writing file (Triangle node file to be overwritten)
  ofstream file;

  // write on .node file
  dummystring = fname->at(0)+".node";
  dummystring = "triangleAfterConvert.node";
  file.open(dummystring.c_str());

  file << (*npoin) << " " << (*ndim) << " " << (*ndof) << " 1\n";
  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   file << IPOIN+1 << " ";
   for(unsigned IA=0; IA<(*ndim); IA++) { file << (*c_XY)(IA,IPOIN) << " "; }
   for(unsigned IA=0; IA<(*ndof); IA++) { file << (*c_Zroe)(IA,IPOIN) << " "; }
   file << nodcod->at(IPOIN) << "\n";
  }

  file.close();

  // write on .poly file
  dummystring = fname->at(0)+".poly";
  dummystring = "triangleAfterConvert.poly";
  file.open(dummystring.c_str());

  file << "0 " << (*ndim) << " 0" << " 1\n";
  file << (*nbfac) << " 1\n";
  for(unsigned IFACE=0; IFACE<(*nbfac); IFACE++) {
   NBND = namebnd.at((*bndfac)(2,IFACE)-1); // c++ indeces start from 0
   if(NBND=="InnerSup" || NBND=="InnerSub") { NBND=10; }
   file << IFACE+1 << " ";
   file << (*bndfac)(0,IFACE) << " " << (*bndfac)(1,IFACE) << " " << NBND << "\n";
  }

  file << "0\n"; // write number of holes   
  file.close();
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setAddress()
{
  unsigned start;
  start = 0;
  c_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
  c_Zroe = new Array2D <double> ((*ndof),(*npoin),&zroe->at(start));
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nelem = MeshData::getInstance().getData <unsigned> ("NELEM");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  nbpoin = MeshData::getInstance().getData <unsigned> ("NBPOIN");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nHole = MeshData::getInstance().getData <unsigned> ("NHOLE");
  zroe = MeshData::getInstance().getData< std::vector<double> >("ZROE");
  coor = MeshData::getInstance().getData< std::vector<double> >("COOR");
  nodcod = MeshData::getInstance().getData< std::vector<int> >("NODCOD");
  celnod = MeshData::getInstance().getData< Array2D<int> >("CELNOD");
  celcel = MeshData::getInstance().getData< Array2D<int> >("CELCEL");
  bndfac = MeshData::getInstance().getData< Array2D<int> >("BNDFAC");
  fname = MeshData::getInstance().getData< std::vector<string> >("FNAME");
}

//----------------------------------------------------------------------------//

void CFmesh2Triangle::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  neshmax = PhysicsData::getInstance().getData <unsigned> ("NESHMAX");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

