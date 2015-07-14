// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/CFmesh2StartingTriangleFreez.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/CrossProd.hh"
#include "MathTools/MinMax.hh"
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
ObjectProvider<CFmesh2StartingTriangleFreez, Converter>
CFmesh2StartingTriangleFreezProv("CFmesh2StartingTriangleFreez");

//----------------------------------------------------------------------------//

CFmesh2StartingTriangleFreez::CFmesh2StartingTriangleFreez(const std::string& objectName) :
  Converter(objectName)
{
  m_meshInputfile = "DummyCFMeshFile";
  addOption("InputFile",&m_meshInputfile,
            "CFmesh file containing the captured solution");
  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

CFmesh2StartingTriangleFreez::~CFmesh2StartingTriangleFreez()
{
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::setup()
{
  LogToScreen(VERBOSE, "CFmesh2StartingTriangleFreez::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "CFmesh2StartingTriangleFreez::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::unsetup()
{
  LogToScreen(VERBOSE, "CFmesh2StartingTriangleFreez::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::configure(OptionMap& cmap, const std::string& prefix)
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

void CFmesh2StartingTriangleFreez::convert()
{
  LogToScreen (INFO, "CFmesh2StartingTriangleFreez::convert()\n");

  logfile.Open(getClassName());

  // read CFmesh format file
  LogToScreen(DEBUG_MIN, "CFmesh2StartingTriangleFreez::reading CFmesh format\n");
  readCFmeshFmt();

  computeEdges();

  computeEdgesFlag();

  // make the transformation from primitive variables to
  // Roe parameter vector variables
  vector <double> m_prim(ndof);
  vector <double> m_zroe(ndof);
  vector <double> m_XY(ndim);
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   for(unsigned i=0; i<ndof; i++) { m_prim.at(i) = prim.at(IPOIN*ndof+i); }
   for(unsigned i=0; i<ndim; i++) { m_XY.at(i) = XY.at(IPOIN*ndim+i); }
   m_prim2param.ptr()->transform(m_prim,m_XY,m_zroe);
   for(unsigned i=0; i<ndof; i++) { zroe.at(IPOIN*ndof+i) = m_zroe.at(i); }
  }

  // write Triangle format files
  LogToScreen(DEBUG_MIN, "CFmesh2StartingTriangleFreez::writing Triangle format\n");
  writeTriangleFmt();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::readCFmeshFmt()
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

  file.open(string(m_meshInputfile).c_str());

  // read number of the first dummy strings
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open(string(m_meshInputfile).c_str());
  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  while(ISKIP<(LSKIP)) { file >> skipver >> dummy;
                          ++ISKIP;      }
  file >> dummy >> ndim;                           // read !NB_DIM    ndim
  file >> dummy >> ndof;                           // read !NB_EQ     ndof
  file >> dummy >> npoin >> dummy;                 // read !NB_NODES  npoin 0
  file >> dummy >> dummy >> dummy;                 // read !NB_STATES nbstates 0
  file >> dummy >> nelem;                          // read !NB_ELEM   nelem
  file >> dummy >> dummy;                          // read !NB_ELEM_TYPES nbelemTypes

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
  while(I<nelem) {
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

  nbfac = 0;
  for(unsigned IFACE=0; IFACE<NCLR; IFACE++) {
   nbfac = nFacB.at(IFACE) + nbfac;
  }
 
  nvt = ndim + 1;
  nholes = 0;

  // resize vectors with the dimensions values read on CFmesh file
  resizeVectors();

  // read mesh points coordinates, mesh points status, the mesh connectivity,
  // the boundary structure and the solution from CFmesh file
  file.open(string(m_meshInputfile).c_str());

  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSIONE cfmeshFmtversion
  LSKIP = 0;
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open(string(m_meshInputfile).c_str());
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
  for(unsigned IELEM=0; IELEM<nelem; IELEM++) {
   for(unsigned IVERT=0; IVERT<nvt; IVERT++) {
    file >> value;
    celnod.at(IELEM*nvt+IVERT) = value + 1;
   }
   for(unsigned IVERT=0; IVERT<nvt; IVERT++) {
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
    bndfac.at(IFACE*3+0) = bndfac0 + 1;
    bndfac.at(IFACE*3+1) = bndfac1 + 1;
    bndfac.at(IFACE*3+2) = IC+1; // c++ indeces start from 0
    ++IFACE;
   }
  }

  // read !LIST_NODE and other dummy strings before the nodal coordinates
  file >> dummy; // read !LIST_NODES

  if(dummy!="!LIST_NODE") { file >> dummy; }

  // read the nodal coordinates
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   for(unsigned IV=0; IV<ndim; IV++) { file >> XY.at(IPOIN*ndim+IV); }
  }
  file >> dummy >>  LIST_STATES;

  if      (LIST_STATES==0) {
   for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
    for(unsigned I=0; I<ndof; I++) { prim.at(IPOIN*ndof+I) = 0; }
   }
  }

  else if (LIST_STATES==1) {
   for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
    for(unsigned I=0; I<ndof; I++) { file >> prim.at(IPOIN*ndof+I); }
   }
  }
  else { cout << "CFmesh2StartingTriangle::error => !LIST_STATES must be 0 or 1 ";
         cout << "in the .CFmesh file\n"; }

  // fill nodcod vector
  setNodcod();

  file.close();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::setNodcod()
{
  unsigned IPOIN;
  int IELEM, IVERT, IBC, KVERT;

  /// create Jcycl object
  Jcycl J;

  // flag with 0 those meshpoints that belong to an element
  for(unsigned IELEM=0; IELEM<nelem; IELEM++) {
   for(unsigned IV=0; IV<nvt; IV++) {
    IPOIN = celnod.at(IELEM*nvt+IV);
    nodcod.at(IPOIN-1) = 0; // c++ indeces start from 0
   }
  }

  // flag with smthg > 0 those meshpoints that belong to bndry face
  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   IELEM = bndfac.at(IFACE*3+0);
   IVERT = bndfac.at(IFACE*3+1);
   IBC = bndfac.at(IFACE*3+2);
   for(unsigned JVERT=0; JVERT<nvt-1; JVERT++) {
    KVERT = J.callJcycl(IVERT+JVERT+1)-1;
    IPOIN = celnod.at((IELEM-1)*nvt+KVERT); // c++ indeces start from 0
    nodcod.at(IPOIN-1) = nodcod.at(IPOIN-1)+1;
   }
  }

  // Count the boundary and hanging nodes
  nhnode = 0;
  nbpoin = 0;
  for(IPOIN=0; IPOIN<npoin; IPOIN++) {
   if      (nodcod.at(IPOIN) == -1) { ++nhnode; }
   else if (nodcod.at(IPOIN) > 0)   { ++nbpoin; }
  }
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::resizeVectors()
{
  prim.resize(ndof*npoin);
  zroe.resize(ndof*npoin);
  XY.resize(ndim*npoin);
  bndfac.resize(3*nbfac);
  nodcod.resize(npoin);
  celnod.resize(3*nelem);
  celcel.resize(3*nelem);
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::computeEdges()
{
  // call Jcycl object
  Jcycl J;

  // call the object computing the maximum and minimum value
  MinMax <double> m;

  // call the objects computing the cross product
  CrossProd <double> c1;
  CrossProd <double> c2;

  int I1, I2;
  int J0, J1, J2;
  int N0, N1, N2, N3;
  int M1, M2;
  int LL;

  vector <double> D0(3);
  vector <double> D1(3);
  vector <double> D2(3);
  vector <double> AA(3);
  vector <double> BB(3);

  nedges = 3 * npoin - nbfac - 3 + 3*nholes;

  // resize celedg and edgnod arrays
  celedg.resize(nvt*nelem);
  edgnod.resize(4*nedges);

  // initialization of celcel pointer
  for(unsigned IE=0; IE<nelem; IE++) {
   for(unsigned IV=0; IV<nvt; IV++) {
    I1 = celnod.at(IE*nvt+(J.callJcycl(IV+2)-1)); // celnod(J.callJcycl(IV+2)-1,IE)
    I2 = celnod.at(IE*nvt+(J.callJcycl(IV+3)-1)); // celnod(J.callJcycl(IV+3)-1,IE)

    celcel.at(IE*nvt+IV) = 0; // celcel(IV,IE)

    for(unsigned IE2=0; IE2<nelem; IE2++) {
     if(IE2!=IE) {
      J0 = celnod.at(IE2*nvt+0); // celnod(0,IE2)
      J1 = celnod.at(IE2*nvt+1); // celnod(1,IE2)
      J2 = celnod.at(IE2*nvt+2); // celnod(2,IE2)
      if( (I1==J0 || I1==J1 || I1==J2) &&
          (I2==J0 || I2==J1 || I2==J2) ) {
       celcel.at(IE*nvt+IV) = IE2+1; // celcel(IV,IE)
      }
     }
     if(celcel.at(IE*nvt+IV)==0) { celcel.at(IE*nvt+IV) = -1; } // celcel(IV,IE)
    }
   }
  }

  unsigned IEDGE=0;

  for(unsigned IELEM=0; IELEM<nelem; IELEM++) {

   // find the node that make up an edge
   unsigned IVERT=0;
   ivert:
    while(IVERT<nvt) {
     int JELEM = celcel.at(IELEM*nvt+IVERT); // celcel(IVERT,IELEM)

     // are we going to add a new edge?
     if(JELEM > IELEM+1 || JELEM<1 || JELEM>nelem) {
      N0 = celnod.at(IELEM*nvt+IVERT); // celnod(IVERT, IELEM)
      // celnod(J.callJcycl(IVERT+2)-1,IELEM)
      N1 = celnod.at(IELEM*nvt+(J.callJcycl(IVERT+2)-1));
      // celnod(J.callJcycl(IVERT+3)-1,IELEM) 
      N2 = celnod.at(IELEM*nvt+(J.callJcycl(IVERT+3)-1)); 
      ++IEDGE;
      celedg.at(IELEM*nvt+IVERT) = IEDGE; // celedg(IVERT,IELEM)

      // boundary edges
      edgnod.at((IEDGE-1)*4+0) = N1; // edgnod(0,IEDGE-1)
      edgnod.at((IEDGE-1)*4+1) = N2; // edgnod(1,IEDGE-1) 

      if(JELEM > IELEM+1 || JELEM <1) {
       edgnod.at((IEDGE-1)*4+2) = IELEM+1; // edgnod(2,IEDGE-1)
       edgnod.at((IEDGE-1)*4+3) = JELEM; // edgnod(3,IEDGE-1) 
      }
      // interior edges
      else {
       // pick up the vertex of cell JELEM that faces the current edge
       for(unsigned JV=0; JV<nvt; JV++) {
        N3 = celnod.at((JELEM-1)*nvt+JV); // celnod(JV,JELEM-1)
        // celnod(J.callJcycl(IVERT+2)-1,JELEM-1) 
        M1 = celnod.at((JELEM-1)*nvt+(J.callJcycl(IVERT+2)-1));
        // celnod(J.callJcycl(IVERT+3)-1,JELEM-1) 
        M2 = celnod.at((JELEM-1)*nvt+(J.callJcycl(IVERT+3)-1)); 
        if(m.min(M1,M2)!=m.min(N1,N2) && m.max(M1,M2)!=m.max(N1,N2)) {
         cout << "CFmesh2StartingTriangleFreez::error => something wrong occurred\n";
         exit(1);
        }
       }
       // let us check on which side is N3
       D0.at(0) = XY.at((N2-1)*PhysicsInfo::getnbDim()+0) - 
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+0); // XY(0,N2-1)-XY(0,N1-1)
       D0.at(1) = XY.at((N2-1)*PhysicsInfo::getnbDim()+1) - 
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+1); // XY(1,N2-1)-XY(1,N1-1)
       D0.at(2) = 0.0;
       D1.at(0) = XY.at((N0-1)*PhysicsInfo::getnbDim()+0) - 
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+0); // XY(0,N0-1)-XY(0,N1-1)
       D1.at(1) = XY.at((N0-1)*PhysicsInfo::getnbDim()+1) -
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+1); // XY(1,N0-1)-XY(1,N1-1)
       D1.at(2) = 0.0;
       D2.at(0) = XY.at((N3-1)*PhysicsInfo::getnbDim()+0) -
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+0); // XY(0,N3-1)-XY(0,N1-1)
       D2.at(1) = XY.at((N3-1)*PhysicsInfo::getnbDim()+1) -
                  XY.at((N1-1)*PhysicsInfo::getnbDim()+1); // XY(1,N3-1)-XY(1,N1-1)
       D2.at(2) = 0.0;

       c1.computeCrossProd(D0,D1);
       AA = c1.getCrossProd();
       c2.computeCrossProd(D0,D2);
       BB = c2.getCrossProd();
       if( AA.at(2) < 0.0 && BB.at(2) > 0.0) {
        // element JELEM is on the left
        edgnod.at((IEDGE-1)*4+2) = JELEM; // edgnod(2,IEDGE-1)
        edgnod.at((IEDGE-1)*4+3) = IELEM+1; // edgnod(3,IEDGE-1)
       }
       else if ( AA.at(2) > 0.0 && BB.at(2)<0.0) {
        // element JELEM is on the right
        edgnod.at((IEDGE-1)*4+2) = IELEM+1; // edgnod(2,IEDGE-1)
        edgnod.at((IEDGE-1)*4+3) = JELEM; // edgnod(3,IEDGE-1)
       }
       else {
        cout << "CFmesh2StartingTriangleFreez::error => something wrong occured in locating elements\n";
       }
      }
     }
     else {
      unsigned JV=0;
      jv:
       while(JV<nvt) {
        LL = celnod.at((JELEM-1)*nvt+JV); // celnod(JV,JELEM-1) 

        for(unsigned IV=0; IV<nvt; IV++) {
          if (celnod.at(IELEM*nvt+IV) == LL) { JV++; goto jv; break; }
        }
        // celedg(IVERT,IELEM) = celedg(JV,JELEM-1)
        celedg.at(IELEM*nvt+IVERT) = celedg.at((JELEM-1)*nvt+JV);
         ++IVERT; goto ivert; break;
       }
       cout << "CFmesh2StartingTriangleFreez::error => Something has gone wrong in cell ";
       cout << IELEM+1 << endl;
       for(unsigned K=0; K<nvt; K++) 
        { cout << celnod.at(IELEM*nvt+K) << " "; } // celnod(K,IELEM)
       cout << endl;
       IVERT++; goto ivert; break;
     }
    ++IVERT;
    }
   }

   if(IEDGE!=nedges) {
    logfile("\nCFmesh2StartingTriangleFreez::warning => problem with the number of edges\n");
    logfile("                                 Should be ",nedges) ;
    logfile(" is ",IEDGE,"\n\n");
    if (nedges < IEDGE) {
     cout << "CFmesh2StartingTriangleFreez::error => a memory allocation problem has occurred\n";
     exit(1);
    }
    else {
     nedges = IEDGE;
    }
   }
   else {
    logfile("CFmesh2StartingTriangleFreez::info => check on edges is all right, found: ");
    logfile(IEDGE,"\n");
   }
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::computeEdgesFlag()
{
  flag.resize(nedges);

  // call objct computing minumum and maximum values
  MinMax <double> m;

  vector <double> D1(2);

  int N0 = 0;
  int I1, I2, N1, N2, J1, J2;
  double D0;

  unsigned IEDGE=0;
  iedge:
   while(IEDGE<nedges) {
    int JELEM = edgnod.at(IEDGE*4+3); // edgnod(3,IEDGE)
    // test of boundary edges
    if(JELEM>nelem || JELEM <1) {
     I1 = edgnod.at(IEDGE*4+0); // edgnod(0,IEDGE)
     I2 = edgnod.at(IEDGE*4+1); // edgnod(1,IEDGE)

     for(unsigned IBFAC=0; IBFAC<nbfac; IBFAC++) {
      J1 = bndfac.at(IBFAC*3+0); // bndfac(0,IBFAC)
      J2 = bndfac.at(IBFAC*3+1); // bndfac(1,IBFAC)

      if(m.min(I1,I2)==m.min(J1,J2) && m.max(I1,I2)==m.max(J1,J2)) {
       flag.at(IEDGE) = bndfac.at(IBFAC*3+2); // bndfac(2,IBFAC)
       ++N0;

       ++ IEDGE;
       goto iedge; break;
      }
     }
     cout << "CFmesh2StartingTriangleFreez::error => cannot match boundary edge " << IEDGE+1;
     cout << "                               \n edge pointer is: \n";
     for(unsigned K=0; K<4; K++) {
      cout << edgnod.at(IEDGE*4+K) << " "; // edgnod(K,IEDGE)
     }
     cout << endl;
     exit(1);
    }
    else {
     N1 = edgnod.at(IEDGE*4+0)-1; // edgnod(0,IEDGE)-1 
     N2 = edgnod.at(IEDGE*4+1)-1; // edgnod(1,IEDGE)-1 
     // sqrt(pow((*XY)(0,N1),2) + pow((*XY)(1,N1),2))
     D1.at(0) = sqrt(pow(XY.at(N1*PhysicsInfo::getnbDim()+0),2) +
                     pow(XY.at(N1*PhysicsInfo::getnbDim()+1),2));
     // sqrt(pow((*XY)(0,N2),2) + pow((*XY)(1,N2),2))
     D1.at(1) = sqrt(pow(XY.at(N2*PhysicsInfo::getnbDim()+0),2) +
                     pow(XY.at(N2*PhysicsInfo::getnbDim()+1),2));
     D0 = m.max(D1.at(0),D1.at(1));
     if(D0 <= 1.02e0) {
      ++N0;
      flag.at(IEDGE) = 999;
     }
    }
    ++IEDGE; goto iedge;
   }

   logfile("CFmesh2StartingTriangleFreez::info => ",N0);
   logfile(" edges have been flagged out of ",nedges,"\n");
   logfile("                               boundary edges are ");
   logfile(nbfac,"\n");
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangleFreez::writeTriangleFmt()
{
  // count the number of boundary nodes
  unsigned inum = 0;

  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  string NBND;

  // writing file (Triangle node file to be overwritten)
  FILE* trianglefile;

  // write on .node file
  dummystring = "na00.node";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%u %s %u",npoin," ",PhysicsInfo::getnbDim());
  fprintf(trianglefile,"%s %u %s"," ",ndof," 1\n");
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   if(nodcod.at(IPOIN)>0) { ++inum; }
   fprintf(trianglefile,"%u %s",IPOIN+1," ");
   for(unsigned IA=0; IA<ndim; IA++)
    { fprintf(trianglefile,"%20.16F %s",XY.at(IPOIN*ndim+IA)," "); }
   for(unsigned IA=0; IA<ndof; IA++)
   { fprintf(trianglefile,"%20.16F %s",zroe.at(IPOIN*ndof+IA)," "); }
   fprintf(trianglefile,"%u %s",nodcod.at(IPOIN),"\n");
  }

  fclose(trianglefile);

  if(inum==0) {
   cout << "CFmesh2StartingTriangleFreez::info => there are not boundary nodes flagged\n";
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
  dummystring = "na00.poly";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%s %u %s", "0 ", ndim, " 0 1\n");
  fprintf(trianglefile,"%u %s",nbfac+N0," 1\n");
  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   NBND = namebnd.at(bndfac.at(IFACE*3+2)-1); // c++ indeces start from 0
   if(NBND=="InnerSup" || NBND=="InnerSub") { NBND="10"; }
   fprintf(trianglefile,"%u %s",IFACE+1," ");
   fprintf(trianglefile,"%i %s %i",bndfac.at(IFACE*3+0)," ",bndfac.at(IFACE*3+1));
   fprintf(trianglefile,"%s %s %s"," ",NBND.c_str()," \n");
  }

  /// freezed edges
  unsigned IEDG=0;
  for(unsigned IFACE=0; IFACE<nedges; IFACE++) {
   unsigned N1 = edgnod.at(IFACE*4+0); // edgnod(0,IFACE)
   unsigned N2 = edgnod.at(IFACE*4+1); // edgnod(1,IFACE)
   if(flag.at(IFACE)==999) {
    ++IEDG;
    fprintf(trianglefile,"%u %s %u %s",nbfac+IEDG," ",N1," ");
    fprintf(trianglefile,"%u %s %i %s",N2," ",flag.at(IFACE),"\n");
   }
  }

  fprintf(trianglefile,"%s","0\n"); // write number of holes   
  fclose(trianglefile);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

