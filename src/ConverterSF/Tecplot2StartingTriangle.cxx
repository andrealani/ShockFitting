// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/Tecplot2StartingTriangle.hh"
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
ObjectProvider<Tecplot2StartingTriangle, Converter>
Tecplot2StartingTriangleProv("Tecplot2StartingTriangle");

//----------------------------------------------------------------------------//

Tecplot2StartingTriangle::Tecplot2StartingTriangle(const std::string& objectName) :
  Converter(objectName)
{
  // first the *.plt and then the *surf.plt
  m_meshInputfile = vector<string>(); 
  addOption("InputFile",&m_meshInputfile,
            "Tecplot files containing the captured solution");
  // it must be set equal to 0 per Perfect Gas model
  m_nbSpecies = 0;
  addOption("nbSpecies",&m_nbSpecies,
            "Specifies the number of chemical species");
  m_tecplotExtraValues = false;
  addOption("extraValuesPrinted",&m_tecplotExtraValues,
            "Specifies if extra values are printed in the tecplot file");

  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

Tecplot2StartingTriangle::~Tecplot2StartingTriangle()
{
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTriangle::setup()
{
  LogToScreen(VERBOSE, "Tecplot2StartingTriangle::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "Tecplot2StartingTriangle::setup() => end\n");
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTriangle::unsetup()
{
  LogToScreen(VERBOSE, "Tecplot2StartingTriangle::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTriangle::configure(OptionMap& cmap, const std::string& prefix)
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

void Tecplot2StartingTriangle::convert()
{
  LogToScreen (INFO, "Tecplot2StartingTriangle::convert()\n");

  // read Tecplot format file
  LogToScreen(DEBUG_MIN, "Tecplot2StartingTriangle::reading Tecplot format\n");
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
  LogToScreen(DEBUG_MIN, "Tecplot2StartingTriangle::writing Triangle format\n");
  writeTriangleFmt();
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTriangle::readTecplotFmt()
{
  // dummy variables
  string dumstring; double dumvalue;
  string::size_type i;
  vector <unsigned> m_bndfac(2);
  unsigned countbnd;

  // number of point for each coloured face
  unsigned npoinFac;

  // map array used to match tecplot id-boundary edges
  // and id-boundary edges
  // it will be resized after
  Array2D <unsigned> array_bnd(1,1);

  // coordinates of the boundary points
  vector <double> XY_bnd(2);

  // reading file
  ifstream file;

  ndof = m_nbSpecies+4;

  file.open(string(m_meshInputfile.at(0)).c_str());

  // read TITLE = ...
  getline(file,dumstring);
  // read VARIABLES = x x x ...
  getline(file,dumstring);
 
  // read ZONE   T="P0 ZONE0 Triag",
  file >> dumstring >> dumstring >> dumstring >> dumstring;

  // read N=....
  file >> dumstring;

  // take the number of points from the string
  i = dumstring.find(string("N="));
  if (i != std::string::npos) { dumstring.erase(i,string("N=").length()); }
  i = dumstring.find(string(","));
  if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

  npoin = atoi(dumstring.c_str());

  // read E=....
  file >> dumstring;
  
  // take the number of elements from the string
  i = dumstring.find(string("E="));
  if (i != std::string::npos) { dumstring.erase(i,string("E=").length()); }
   i = dumstring.find(string(","));
  if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

  nelem = atoi(dumstring.c_str());

  // read F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0
  getline(file,dumstring);

  prim.resize(ndof*npoin);
  zroe.resize(ndof*npoin);
  XY.resize(PhysicsInfo::getnbDim()*npoin);
  nodcod.resize(npoin);

  for(unsigned IPOIN=0;IPOIN<npoin;IPOIN++) {
   // read the nodal coordinates
   for(unsigned IV=0; IV<PhysicsInfo::getnbDim(); IV++)
    { file >> XY.at(IPOIN*PhysicsInfo::getnbDim()+IV); }
   // read the nodal states
   for(unsigned I=0; I<ndof; I++) { file >> prim.at(IPOIN*ndof+I); }
   if(m_tecplotExtraValues) {
    for(unsigned I=0; I<4; I++) { file >> dumvalue; }
   }
  } 

  // the list of elements will not be read, celnod and celcel will
  // be generated calling triangle
 
  file.close(); // close *.plt

  file.open(string(m_meshInputfile.at(1)).c_str()); // open *surf.plt 

  // read TITLE = ...
  getline(file,dumstring);
  // read VARIABLES = x x x ...
  getline(file,dumstring);

  unsigned IFACE=0;
  nbfac = 0;
  while(!file.eof()) {
   
   // read ZONE
   file >> dumstring;

   if(dumstring=="ZONE" || dumstring=="ZONE " || dumstring=="\nZONE") {
    // read "N=.."
    file >> dumstring;

    i = dumstring.find(string("N="));
    if (i != std::string::npos) { dumstring.erase(i,string("N=").length()); }
    i = dumstring.find(string(","));
    if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    npoinFac = atoi(dumstring.c_str());

    array_bnd.resize(2,npoinFac);

    // read "T=.."
    file >> dumstring;
    i = dumstring.find(string("T=\""));
    if (i != std::string::npos) { dumstring.erase(i,string("T=\"").length()); }
    i = dumstring.find(string(","));
    if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    namebnd.resize(IFACE+1);

    namebnd.at(IFACE)=atoi(dumstring.c_str());

    // read  "TR 0"," 
    file >> dumstring >> dumstring;

    // read "E=.."
    file >> dumstring;

    i = dumstring.find(string("E="));
    if (i != std::string::npos) { dumstring.erase(i,string("E=").length()); }
    i = dumstring.find(string(","));
    if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    nFacB.resize(IFACE+1);
    nFacB.at(IFACE) = atoi(dumstring.c_str());

    bndfac.resize(3*(bndfac.size()+nFacB.at(IFACE)));
 
    // read F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0
    file >> dumstring  >> dumstring  >> dumstring;

    // read boundary nodes values
    // and assign the correct IDnumber to the boundary points
    countbnd=0;
    for(unsigned IBPOIN=0;IBPOIN<npoinFac;IBPOIN++) {
     for(unsigned IV=0;IV<PhysicsInfo::getnbDim();IV++) { 
      file >> XY_bnd.at(IV);
     }
     for(unsigned IPOIN=0;IPOIN<npoin;IPOIN++) {
      if(XY_bnd.at(0)==XY.at(IPOIN*PhysicsInfo::getnbDim()+0) &&
         XY_bnd.at(1)==XY.at(IPOIN*PhysicsInfo::getnbDim()+1)) {
       array_bnd(0,IBPOIN)=IPOIN+1;
       array_bnd(1,IBPOIN)=IBPOIN+1;
       countbnd++;
      }
     }
     for(unsigned IV=0;IV<ndof;IV++) { file >> dumvalue; }
     if(m_tecplotExtraValues) {
      for(unsigned IV=0; IV<4; IV++) { file >> dumvalue; }
     }
    }
    if(countbnd!=npoinFac) {
     cout << "Tecplot2StartingTriangle::error => no match between\n";
     cout << "  coordinates of boundary points\n";
     cout << "countbnd= "<<countbnd<<" # boundary points= "<<npoinFac<<endl;
     exit(1);
    }
 
    // read geometry entities list that is boundary data
    countbnd=0;
    for(unsigned IB=nbfac;IB<nFacB.at(IFACE)+nbfac;IB++) {
     for(int IV=0;IV<2;IV++) { 
      file >> m_bndfac.at(IV);
      for(unsigned IBPOIN=0;IBPOIN<npoinFac;IBPOIN++) {
       if(m_bndfac.at(IV)==array_bnd(1,IBPOIN)) {
        bndfac.at(IB*3+IV)=array_bnd(0,IBPOIN); countbnd++; break; }
      } // for IBPOIN
     } // for IV=1:ndim
     bndfac.at(IB*3+2) = namebnd.at(IFACE);
    } 
    if(countbnd!=nFacB.at(IFACE)*2) {
     cout << "Tecplot2StartingTriangle::error => wrong map for\n";
     cout << "tecplot ID boundary edges and ID boundary edges\n";
     cout << "countbnd= "<<countbnd<<" # edge points= "<<nFacB.at(IFACE)*2<<endl;
     exit(1);
    }

    // update the total number of boundary edges
    nbfac = nFacB.at(IFACE) + nbfac;

    IFACE++;
   } // if dumstring=ZONE

  } // while !eof

  NCLR = IFACE;
/*nbfac=0;
for(unsigned IF=0;IF<NCLR;IF++) {
for(unsigned IB=nbfac;IB<nFacB.at(IF)+nbfac;IB++) {
cout << "IB " << IB << endl;
for(int IV=0;IV<3;IV++) {
cout << bndfac.at(IB*3+IV) << " ";
}
cout << endl;
}
nbfac = nFacB.at(IFACE) + nbfac;
cout << nbfac << endl;
}*/

  file.close(); // close *surf.plt

  // fill nodcod vector
  setNodcod();
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTriangle::setNodcod()
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

void Tecplot2StartingTriangle::writeTriangleFmt()
{
  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  stringstream NBND;

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

  unsigned m_nbfac=0;
  for(unsigned IFACE=0;IFACE<NCLR;IFACE++)  {
   for(unsigned IB=m_nbfac;IB<nFacB.at(IFACE)+m_nbfac;IB++) {
    NBND.str(string());
    NBND << namebnd.at(IFACE); // c++ indeces start from 0
    if(NBND=="InnerSup" || NBND=="InnerSub") { NBND << "10"; }
    fprintf(trianglefile,"%5u",IB+1);
    fprintf(trianglefile,"%5i %5i",bndfac.at(IB*3+0),bndfac.at(IB*3+1));
    fprintf(trianglefile,"%3s %s",NBND.str().c_str()," \n");
   }
   m_nbfac=m_nbfac+nFacB.at(IFACE);
  }

  fprintf(trianglefile,"%s","0\n"); // write number of holes   
  fclose(trianglefile);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

