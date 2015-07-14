// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/CFmesh2StartingTriangle.hh"
#include "Framework/Log.hh"
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
ObjectProvider<CFmesh2StartingTriangle, Converter>
CFmesh2StartingTriangleProv("CFmesh2StartingTriangle");

//----------------------------------------------------------------------------//

CFmesh2StartingTriangle::CFmesh2StartingTriangle(const std::string& objectName) :
  Converter(objectName)
{
  m_meshInputfile = "";
  addOption("InputFile",&m_meshInputfile,
            "CFmesh file containing the captured solution");
  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

CFmesh2StartingTriangle::~CFmesh2StartingTriangle()
{
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::setup()
{
  LogToScreen(VERBOSE, "CFmesh2StartingTriangle::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "CFmesh2StartingTriangle::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::unsetup()
{
  LogToScreen(VERBOSE, "CFmesh2StartingTriangle::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::configure(OptionMap& cmap, const std::string& prefix)
{
  Converter::configure(cmap, prefix);

  // assign strings on input.case file to variable transformer object
  m_prim2param.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;

  if (ConfigFileReader::isFirstConfig()) {
   m_prim2param.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                             getProvider(m_prim2param.name())
                             ->create(m_prim2param.name()));
  }
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::convert()
{
  LogToScreen (INFO, "CFmesh2StartingTriangle::convert()\n");

  // read CFmesh format file
  LogToScreen(DEBUG_MIN, "CFmesh2StartingTriangle::reading CFmesh format\n");
  readCFmeshFmt();

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
  LogToScreen(DEBUG_MIN, "CFmesh2StartingTriangle::writing Triangle format\n");
  writeTriangleFmt();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::readCFmeshFmt()
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
  // read !CFMESH_FORMAT_VERSION cfmeshFmtversion
  while(ISKIP<(LSKIP)) { file >> skipver >> dummy;
                          ++ISKIP;      }
  file >> dummy >> ndim;                              // read !NB_DIM    ndim
  file >> dummy >> ndof;                              // read !NB_EQ     ndof
  file >> dummy >> npoin >> dummy;                    // read !NB_NODES  npoin 0
  file >> dummy >> dummy >> dummy;                    // read !NB_STATES nbstates 0
  file >> dummy >> nelem;                             // read !NB_ELEM   nelem
  file >> dummy >> dummy;                             // read !NB_ELEM_TYPES nbelemTypes

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
   for(unsigned IV=0; IV<ndim; IV++) { file >> XY.at(IPOIN*ndim+IV);}
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
  else { cout << "CFmesh2Triangle::error => !LIST_STATES must be 0 or 1 ";
         cout << "in the .CFmesh file\n"; }

  // fill nodcod vector
  setNodcod();

  file.close();
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::setNodcod()
{
  int IPOIN, help;

  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   for(unsigned I=0; I<nvt-1; I++) {
    IPOIN = bndfac.at(IFACE*3+I);
    help = nodcod.at(IPOIN-1);
    nodcod.at(IPOIN-1) = help+1;
   }
  }
/*
  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   for(unsigned I=0; I<nvt-1; I++) {
    IPOIN = bndfac.at(IFACE*3+I);
    nodcod.at(IPOIN-1) = bndfac.at(IFACE*3+2);
   }
  }
*/
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::resizeVectors()
{
  prim.resize(ndof*npoin);
  zroe.resize(ndof*npoin);
  XY.resize(ndim*npoin);
  bndfac.resize(3*nbfac);
  nodcod.resize(npoin);
}

//----------------------------------------------------------------------------//

void CFmesh2StartingTriangle::writeTriangleFmt()
{
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
   fprintf(trianglefile,"%u %s",IPOIN+1," ");
   for(unsigned IA=0; IA<ndim; IA++)
    { fprintf(trianglefile,"%20.16F %s",XY.at(IPOIN*ndim+IA)," "); }
   for(unsigned IA=0; IA<ndof; IA++)
   { fprintf(trianglefile,"%20.16F %s",zroe.at(IPOIN*ndof+IA)," "); }
   fprintf(trianglefile,"%u %s",nodcod.at(IPOIN),"\n");
  }

  fclose(trianglefile);
 
  // write on .poly file
  dummystring = "na00.poly";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%s %u %s", "0 ", ndim, " 0 1\n");
  fprintf(trianglefile,"%u %s",nbfac," 1\n");
  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   NBND = namebnd.at(bndfac.at(IFACE*3+2)-1); // c++ indeces start from 0
   fprintf(trianglefile,"%u %s",IFACE+1," ");
   fprintf(trianglefile,"%i %s %i",bndfac.at(IFACE*3+0)," ",bndfac.at(IFACE*3+1));
   fprintf(trianglefile,"%s %s %s"," ",NBND.c_str()," \n");
  }

  fprintf(trianglefile,"%s","0\n"); // write number of holes   
  fclose(trianglefile);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

