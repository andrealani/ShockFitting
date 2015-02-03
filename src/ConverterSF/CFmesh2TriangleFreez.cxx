// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/CFmesh2TriangleFreez.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Jcycl.hh"
#include "MathTools/CrossProd.hh"
#include "MathTools/MinMax.hh"
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
ObjectProvider<CFmesh2TriangleFreez, Converter>
CFmesh2TriangleFreezProv("CFmesh2TriangleFreez");

//----------------------------------------------------------------------------//

CFmesh2TriangleFreez::CFmesh2TriangleFreez(const std::string& objectName) :
  Converter(objectName)
{
  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

CFmesh2TriangleFreez::~CFmesh2TriangleFreez()
{
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::setup()
{
  LogToScreen(VERBOSE, "CFmesh2TriangleFreez::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "CFmesh2TriangleFreez::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::unsetup()
{
  LogToScreen(VERBOSE, "CFmesh2TriangleFreez::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::configure(OptionMap& cmap, const std::string& prefix)
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

void CFmesh2TriangleFreez::convert()
{
  LogToScreen (INFO, "CFmesh2TriangleFreez::convert()\n");

  setMeshData();

  setPhysicsData();

  logfile.Open(getClassName());

  // read CFmesh format file
  LogToScreen(DEBUG_MIN, "CFmesh2TriangleFreez::reading CFmesh format\n");
  readCFmeshFmt();

  // read nodcod vector from .node file
  readNodCod();

  // compute edges
  computeEdges();

  // compute edges flag
  computeEdgesFlag();

  // make the transformation from primitive variables to 
  // Roe parameter vector variables
  m_prim2param.ptr()->transform(); 

  if (MeshData::getInstance().getVersion()=="original")
  {
   // write Triangle format files
   LogToScreen(DEBUG_MIN, "CFmesh2TriangleFreez::writing Triangle format\n");
   writeTriangleFmt(); 
  }

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::readCFmeshFmt()
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
  else { cout << "CFmesh2TriangleFreez::error => !LIST_STATES must be 0 or 1 ";
         cout << "in cfout.CFmesh file\n"; }

  file.close();
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::resizeVectors()
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

void CFmesh2TriangleFreez::readNodCod()
{
  // string for .node file
  string dummystring;

  // dummy variable
  double dum;

  // fstream variable reading .node file
  ifstream trianglefile;

  dummystring = fname->str()+".1.node";

  // open .node file
  trianglefile.open(dummystring.c_str());

  unsigned startNodcod = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                            PhysicsInfo::getnbShPointsMax();

  trianglefile >> dum >> dum >> dum >> dum;
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   trianglefile >> dum;
   for(unsigned IV=0; IV<PhysicsInfo::getnbDim(); IV++) {
    trianglefile >> dum;
   }
   for(unsigned IV=0; IV<(*ndof); IV++) {
    trianglefile >> dum;
   }
   trianglefile >> nodcod->at(startNodcod+IPOIN-1);
  }
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::computeEdges()
{

  // call Jcycl object 
  Jcycl J;

  // call the object computing the maximum and minimum value
  MinMax <double> m;

  // call the objects computing the cross product
  CrossProd <double> c1;
  CrossProd <double> c2;
  
  unsigned I1, I2;
  unsigned J0, J1, J2;
  unsigned N0, N1, N2, N3;
  unsigned M1, M2;
  int LL;

  vector <double> D0(3);
  vector <double> D1(3);
  vector <double> D2(3);
  vector <double> AA(3);
  vector <double> BB(3);

  nedges = 3 * npoin->at(1) - nbfac->at(1) - 3 + 3*nhole->at(1);

  // resize celedg and edgnod arrays
  celedg = new Array2D <int> ((*nvt),nelem->at(1));
  edgnod = new Array2D <int> (4,nedges);

  // initialization of celcel pointer
  for(unsigned IE=0; IE<nelem->at(1); IE++) {
   for(unsigned IV=0; IV<(*nvt); IV++) {
    I1 = (*celnod)(J.callJcycl(IV+2)-1,IE);
    I2 = (*celnod)(J.callJcycl(IV+3)-1,IE);
    (*celcel)(IV,IE) = 0;

    for(unsigned IE2=0; IE2<nelem->at(1); IE2++) {
     if(IE2!=IE) {
      J0 = (*celnod)(0,IE2);
      J1 = (*celnod)(1,IE2);
      J2 = (*celnod)(2,IE2);
      if( (I1==J0 || I1==J1 || I1==J2) &&
          (I2==J0 || I2==J1 || I2==J2) ) {
       (*celcel)(IV,IE) = IE2+1;
      }
     }
     if((*celcel)(IV,IE)==0) { (*celcel)(IV,IE) = -1; }
    }
   }
  }

  unsigned IEDGE=0;
 
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {

   // find the node that make up an edge
   unsigned IVERT=0;
   ivert:
    while(IVERT<(*nvt)) {
     unsigned JELEM = (*celcel)(IVERT,IELEM);

     // are we going to add a new edge?
     if(JELEM > IELEM+1 || JELEM<1 || JELEM>nelem->at(1)) {
      N0 = (*celnod)(IVERT, IELEM);
      N1 = (*celnod)(J.callJcycl(IVERT+2)-1,IELEM);
      N2 = (*celnod)(J.callJcycl(IVERT+3)-1,IELEM);
      ++IEDGE;
      (*celedg)(IVERT,IELEM) = IEDGE;

      // boundary edges
      (*edgnod)(0,IEDGE-1) = N1; // c++ indeces start from 0
      (*edgnod)(1,IEDGE-1) = N2; // c++ indeces start from 0

      if(JELEM > IELEM+1 || JELEM <1) {
       (*edgnod)(2,IEDGE-1) = IELEM+1; // c++ indeces start from 0
       (*edgnod)(3,IEDGE-1) = JELEM; // c++ indeces start from 0
      }
      // interior edges
      else {
       // pick up the vertex of cell JELEM that faces the current edge
       for(unsigned JV=0; JV<(*nvt); JV++) {
        N3 = (*celnod)(JV,JELEM-1); // c++ indeces start from 0
        M1 = (*celnod)(J.callJcycl(IVERT+2)-1,JELEM-1); // c++ indeces start from 0
        M2 = (*celnod)(J.callJcycl(IVERT+3)-1,JELEM-1); // c++ indeces start from 0
        if(m.min(M1,M2)!=m.min(N1,N2) && m.max(M1,M2)!=m.max(N1,N2)) {
         cout << "CFmesh2TriangleFreez::error => something wrong occurred\n";
         exit(1);
        } 
       }
       // let us check on which side is N3
       D0.at(0) = (*XY)(0,N2-1)-(*XY)(0,N1-1);
       D0.at(1) = (*XY)(1,N2-1)-(*XY)(1,N1-1);
       D0.at(2) = 0.0;
       D1.at(0) = (*XY)(0,N0-1)-(*XY)(0,N1-1);
       D1.at(1) = (*XY)(1,N0-1)-(*XY)(1,N1-1);
       D1.at(2) = 0.0;
       D2.at(0) = (*XY)(0,N3-1)-(*XY)(0,N1-1);
       D2.at(1) = (*XY)(1,N3-1)-(*XY)(1,N1-1);
       D2.at(2) = 0.0;

       c1.computeCrossProd(D0,D1);
       AA = c1.getCrossProd();
       c2.computeCrossProd(D0,D2);
       BB = c2.getCrossProd();
       if( AA.at(2) < 0.0 && BB.at(2) > 0.0) {
        // element JELEM is on the left
        (*edgnod)(2,IEDGE-1) = JELEM;
        (*edgnod)(3,IEDGE-1) = IELEM+1;
       }
       else if ( AA.at(2) > 0.0 && BB.at(2)<0.0) {
        // element JELEM is on the right
        (*edgnod)(2,IEDGE-1) = IELEM+1;
        (*edgnod)(3,IEDGE-1) = JELEM;
       }
       else { 
        cout << "CFmesh2TriangleFreez::error => something wrong occured in locating elements\n";
       }
      }
     }
     else {
      unsigned JV=0;
      jv:
       while(JV<(*nvt)) {
        LL = (*celnod)(JV,JELEM-1); // c++ indeces start from 0
        for(unsigned IV=0; IV<(*nvt); IV++) {
          if ((*celnod)(IV,IELEM) == LL) { JV++; goto jv; break; }
        }
        (*celedg)(IVERT,IELEM) = (*celedg)(JV,JELEM-1); // c++ indeces start from 0
         ++IVERT; goto ivert; break; 
       }
       cout << "CFmesh2TriangleFreez::error => Something has gone wrong in cell ";
       cout << IELEM+1 << endl;
       for(unsigned K=0; K<(*nvt); K++) { cout << (*celnod)(K,IELEM) << " "; }
       cout << endl;
       IVERT++; goto ivert; break; 
     }
    ++IVERT;
    }
   }

   if(IEDGE!=nedges) {
    logfile("\nCFmesh2TriangleFreez::warning => problem with the number of edges\n");
    logfile("                                 Should be ",nedges) ;
    logfile(" is ",IEDGE,"\n\n"); 
    if (nedges < IEDGE) {
     cout << "CFmesh2TriangleFreez::error => a memory allocation problem has occurred\n";
     exit(1);
    }
    else {
     nedges = IEDGE;
    }
   }
   else { 
    logfile("CFmesh2TriangleFreez::info => check on edges is all right, found: ");
    logfile(IEDGE,"\n");
   }
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::computeEdgesFlag()
{
  flag.resize(nedges);

  // call objct computing minumum and maximum values
  MinMax <double> m;

  vector <double> D1(2);
  vector <double> XYC(2);

  XYC.at(0) = 1.0;
  XYC.at(1) = 0.0;

  unsigned N0 = 0;
  unsigned I1, I2, N1, N2, J1, J2;
  double D0;


  unsigned IEDGE=0;
  iedge:
   while(IEDGE<nedges) {
    unsigned JELEM = (*edgnod)(3,IEDGE);

    // test of boundary edges
    if(JELEM>nelem->at(1) || JELEM <1) {
     I1 = (*edgnod)(0,IEDGE);
     I2 = (*edgnod)(1,IEDGE);

     for(unsigned IBFAC=0; IBFAC<nbfac->at(1); IBFAC++) {
      J1 = (*bndfac)(0,IBFAC);
      J2 = (*bndfac)(1,IBFAC);
      if(m.min(I1,I2)==m.min(J1,J2) && m.max(I1,I2)==m.max(J1,J2)) {
       flag.at(IEDGE) = (*bndfac)(2,IBFAC);
       ++N0;
       ++ IEDGE;
       goto iedge; break;
      }
     }
     cout << "CFmesh2TriangleFreez::error => cannot match boundary edge " << IEDGE+1;
     cout << "                               \n edge pointer is: \n";
     for(unsigned K=0; K<4; K++) {
      cout << (*edgnod)(K,IEDGE) << " "; }
     cout << endl;
     exit(1);
    }
    else {
     N1 = (*edgnod)(0,IEDGE)-1; // c++ indeces start from 0
     N2 = (*edgnod)(1,IEDGE)-1; // c++ indeces start from 0
     D1.at(0) = sqrt(pow((*XY)(0,N1)-XYC.at(0),2) + pow((*XY)(1,N1)-XYC.at(1),2));
     D1.at(1) = sqrt(pow((*XY)(0,N2)-XYC.at(0),2) + pow((*XY)(1,N2)-XYC.at(1),2));
     D0 = m.max(D1.at(0),D1.at(1));
     if(D0 <= 1.02e0) {
      ++N0;
      flag.at(IEDGE) = 999; 
     }
    }
    ++IEDGE; goto iedge;
   }

   logfile("CFmesh2TriangleFreez::info => ",N0);
   logfile(" edges have been flagged out of ",nedges,"\n");
   logfile("                               boundary edges are ");
   logfile(nbfac->at(1),"\n");
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::writeTriangleFmt()
{
  // count the number of boundary nodes
  unsigned inum=0;

  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  string NBND;
 
  // writing file (Triangle node file to be overwritten)
  FILE* trianglefile;

  // take new nodcode values from the new nodcod vector of the shocked mesh
  // in the fortran version this new vector is referred to index 1 (NODCOD(1))
  // here it is pushed back to the nodcod of the background mesh
  unsigned startNodcod = 
   npoin->at(0) +2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();

  // write on .node file
  dummystring = fname->str()+".1.node";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%u %s %u",npoin->at(1)," ",PhysicsInfo::getnbDim());
  fprintf(trianglefile,"%s %u %s"," ",(*ndof)," 1\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {

   // count the number of boundary nodes
   if(nodcod->at(startNodcod+IPOIN-1)>0) {
    ++inum;
   }
   fprintf(trianglefile,"%u %s",IPOIN+1," ");
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++)
    { fprintf(trianglefile,"%.16F %s",(*XY)(IA,IPOIN)," "); }
   for(unsigned IA=0; IA<(*ndof); IA++) 
   { fprintf(trianglefile,"%.16F %s",(*zroe)(IA,IPOIN)," "); }
   fprintf(trianglefile,"%u %s",nodcod->at(startNodcod+IPOIN-1),"\n");
  }

  fclose(trianglefile);

  if(inum==0) {
   cout << "CFmesh2TriangleFreez::info => there are not boundary nodes flagged\n";
   exit(1);
  }
  else {
   logfile(inum," nodes have been flagged\n");
  }

  // compute number of edges
  unsigned N0 =0; 
  for(unsigned IFACE=0; IFACE<nedges; IFACE++) {
   if(flag.at(IFACE)==999) { ++N0; }
  }
  logfile(N0," internal edges have been flagged\n"); 

  // write on .poly file
  dummystring = fname->str()+".1.poly";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%s %u %s", "0 ", PhysicsInfo::getnbDim(), " 0 1\n");
  fprintf(trianglefile,"%u %s",nbfac->at(1)+N0," 1\n");
  for(unsigned IFACE=0; IFACE<nbfac->at(1); IFACE++) {
   NBND = namebnd.at((*bndfac)(2,IFACE)-1); // c++ indeces start from 0
   if(NBND=="InnerSup" || NBND=="InnerSub") { NBND="10"; }
   fprintf(trianglefile,"%u %s",IFACE+1," ");
   fprintf(trianglefile,"%i %s %i",(*bndfac)(0,IFACE)," ",(*bndfac)(1,IFACE));
   fprintf(trianglefile,"%s %s %s"," ",NBND.c_str()," \n");
  }

  // freezed edges
  unsigned IEDG=0;
  for(unsigned IFACE=0; IFACE<nedges; IFACE++) {
   unsigned N1 = (*edgnod)(0,IFACE);
   unsigned N2 = (*edgnod)(1,IFACE);
   if(flag.at(IFACE)==999) { 
    ++IEDG;
    fprintf(trianglefile,"%u %s %u %s",nbfac->at(1)+IEDG," ",N1," ");
    fprintf(trianglefile,"%u %s %i %s",N2," ",flag.at(IFACE),"\n");
   }
  }

  fprintf(trianglefile,"%s","0\n"); // write number of holes   
  fclose(trianglefile);
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::freeArray()
{
  delete XY; delete zroe; delete celnod;
  delete celcel; delete bndfac;
  delete edgnod; delete celedg;
}

//----------------------------------------------------------------------------//

void CFmesh2TriangleFreez::setMeshData()
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

void CFmesh2TriangleFreez::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

