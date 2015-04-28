// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/TecplotFVM2StartingTriangle.hh"
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
ObjectProvider<TecplotFVM2StartingTriangle, Converter>
TecplotFVM2StartingTriangleProv("TecplotFVM2StartingTriangle");

//----------------------------------------------------------------------------//

TecplotFVM2StartingTriangle::TecplotFVM2StartingTriangle(const std::string& objectName) :
  Converter(objectName)
{
  // first the *.CFmesh and the *plt files
  m_meshInputfile = vector<string>(); 
  addOption("InputFile",&m_meshInputfile,
            "Tecplot files containing the captured solution");
  // it must be set equal to 0 per Perfect Gas model
  m_nbSpecies = 0;
  addOption("nbSpecies",&m_nbSpecies,
            "Specifies the number of chemical species");
  m_tecplotExtraValues = true;
  addOption("extraValuesPrinted",&m_tecplotExtraValues,
            "Specifies if extra values are printed in the tecplot file");

  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

TecplotFVM2StartingTriangle::~TecplotFVM2StartingTriangle()
{
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::setup()
{
  LogToScreen(VERBOSE, "TecplotFVM2StartingTriangle::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "TecplotFVM2StartingTriangle::setup() => end\n");
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::unsetup()
{
  LogToScreen(VERBOSE, "TecplotFVM2StartingTriangle::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::configure(OptionMap& cmap, const std::string& prefix)
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

void TecplotFVM2StartingTriangle::convert()
{
  LogToScreen (INFO, "TecplotFVM2StartingTriangle::convert()\n");

  // read Tecplot format file
  LogToScreen(DEBUG_MIN, "TecplotFVM2StartingTriangle::reading Tecplot format\n");
  readTecplotFmt();

  // make the transformation from primitive variables to
  // Roe parameter vector variables
  vector <double> m_prim(ndof);
  vector <double> m_zroe(ndof);
  vector <double> m_XY(PhysicsInfo::getnbDim());
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   for(unsigned i=0; i<ndof; i++) 
    { m_prim.at(i) = prim.at(IPOIN*ndof+i); }
   for(unsigned i=0; i<PhysicsInfo::getnbDim(); i++)
    { m_XY.at(i) = XY.at(IPOIN*PhysicsInfo::getnbDim()+i); }
   m_prim2param.ptr()->transform(&m_prim,&m_XY,&m_zroe);
   for(unsigned i=0; i<ndof; i++) { zroe.at(IPOIN*ndof+i) = m_zroe.at(i); }
  }

  // write Triangle format files
  LogToScreen(DEBUG_MIN, "TecplotFVM2StartingTriangle::writing Triangle format\n");
  writeTriangleFmt();
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::readTecplotFmt()
{
  // dummy variables
  string dumstring; double dumvalue;
  string::size_type i;
  vector <unsigned> m_bndfac(2);
  unsigned LSKIP=0; unsigned ISKIP=0; unsigned value;
  string skipver = "dummy"; string dummy; string strvalue;

  // boundary colours
  unsigned NCLR;

  // map array used to match tecplot id-boundary edges
  // and id-boundary edges
  // it will be resized after
  Array2D <unsigned> array_bnd(1,1);

  // coordinates of the boundary points
  vector <double> XY_bnd(2);

  // reading file
  ifstream file;

  // read the connectivity from the CFmesh file
  // while the LIST_STATE will be read from the tecplot file
 
  file.open(string(m_meshInputfile.at(0)).c_str());

  // read number of the first dummy strings
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open(string(m_meshInputfile.at(0)).c_str());

  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSION cfmeshFmtversion
  while(ISKIP<(LSKIP)) { file >> skipver >> dummy;
                          ++ISKIP;      }
  file >> dummy >> ndim;                          // read !NB_DIM    ndim
  file >> dummy >> ndof;                          // read !NB_EQ     ndof
  file >> dummy >> npoin >> dummy;                // read !NB_NODES  npoin 0
  file >> dummy >> dummy >> dummy;                // read !NB_STATES nbstates 0
  file >> dummy >> nelem;                         // read !NB_ELEM   nelem
  file >> dummy >> dummy;                         // read !NB_ELEM_TYPES nbelemTypes

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
   file >> dummy >> dummy >> dummy >> dummy;
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
    for(unsigned i=0; i<5; i++) {file >> dummy; }
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
  file.open(string(m_meshInputfile.at(0)).c_str());

  // read !COOLFLUID_VERSION coolfluidVersion
  // read !COOLFLUID_SVNVERSION coolfluidSvnVersion
  // read !CFMESH_FORMAT_VERSION cfmeshFmtversion
  LSKIP = 0;
  do {
   file >> skipver >> dummy;
   if (skipver == "!NB_DIM") { file.close(); break; }
   ++LSKIP; } while(skipver != "!NB_DIM");

  file.open(string(m_meshInputfile.at(0)).c_str());
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
   file >> value;
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
    file >> dummy >> dummy >> bndfac0 >> bndfac1 >> dummy;
    bndfac.at(IFACE*3+0) = bndfac0 + 1;
    bndfac.at(IFACE*3+1) = bndfac1 + 1;
    bndfac.at(IFACE*3+2) = IC+1; // c++ indeces start from 0
    ++IFACE;
   }
  }

  file.close();


  ndof = m_nbSpecies+4;

  // open the *plt file
  file.open(string(m_meshInputfile.at(1)).c_str());

  // read TITLE = ...
  getline(file,dumstring);
  // read VARIABLES = x x x ...
  getline(file,dumstring);
  // read ZONE   T="P0 ZONE0 Triag",.. 
  getline(file,dumstring);

  // read the nodal coordinates and the list state
  for(unsigned IPOIN=0;IPOIN<npoin;IPOIN++) {
   for(unsigned IV=0; IV<ndim; IV++)
    { file >> XY.at(IPOIN*ndim+IV); }
   for(unsigned I=0; I<ndof; I++) { file >> prim.at(IPOIN*ndof+I); }
   if(m_tecplotExtraValues) {
    for(unsigned I=0; I<4; I++) { file >> dumvalue; }
   }
  } 

  file.close(); // close *.plt

  // fill nodcod vector
  setNodcod();
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::setNodcod()
{
  int IPOIN;

  for(unsigned IB=0; IB<nbfac; IB++) {
   for(unsigned I=0; I<2; I++) {
    IPOIN = bndfac.at(IB*3+I);
    nodcod.at(IPOIN-1) = 2;
   }
  }
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::resizeVectors()
{
  prim.resize(ndof*npoin);
  zroe.resize(ndof*npoin);
  XY.resize(ndim*npoin);
  bndfac.resize(3*nbfac);
  nodcod.resize(npoin);
}

//----------------------------------------------------------------------------//

void TecplotFVM2StartingTriangle::writeTriangleFmt()
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
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++)
    { fprintf(trianglefile,"%20.16F %s",XY.at(IPOIN*PhysicsInfo::getnbDim()+IA)," "); }
   for(unsigned IA=0; IA<ndof; IA++)
   { fprintf(trianglefile,"%20.16F %s",zroe.at(IPOIN*ndof+IA)," "); }
   fprintf(trianglefile,"%u %s",nodcod.at(IPOIN),"\n");
  }

  fclose(trianglefile);
 
  // write on .poly file
  dummystring = "na00.poly";

  trianglefile = fopen(dummystring.c_str(),"w");

  fprintf(trianglefile,"%s %u %s", "0 ", PhysicsInfo::getnbDim(), " 0 1\n");
  fprintf(trianglefile,"%u %s",nbfac," 1\n");
  for(unsigned IFACE=0; IFACE<nbfac; IFACE++) {
   NBND = namebnd.at(bndfac.at(IFACE*3+2)-1); // c++ indeces start from 0
   if(NBND=="InnerSup" || NBND=="InnerSub") { NBND="10"; }
   fprintf(trianglefile,"%u %s",IFACE+1," ");
   fprintf(trianglefile,"%i %s %i",bndfac.at(IFACE*3+0)," ",bndfac.at(IFACE*3+1));
   fprintf(trianglefile,"%s %s %s"," ",NBND.c_str()," \n");
  }

  fprintf(trianglefile,"%s","0\n"); // write number of holes   
  fclose(trianglefile);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

