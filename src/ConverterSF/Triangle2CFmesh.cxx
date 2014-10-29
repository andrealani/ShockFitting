// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/Triangle2CFmesh.hh"
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
ObjectProvider<Triangle2CFmesh, Converter>
triangle2CFmeshProv("Triangle2CFmesh");

//--------------------------------------------------------------------------//

Triangle2CFmesh::Triangle2CFmesh(const std::string& objectName) :
  Converter(objectName)
{
  m_param2prim.name() = "DummyVariableTransformer";
}

//----------------------------------------------------------------------------//

Triangle2CFmesh::~Triangle2CFmesh()
{
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setup()
{
  LogToScreen (VERBOSE, "Triangle2CFmesh::setup() => start\n");

  m_param2prim.ptr()->setup();

  LogToScreen (VERBOSE, "Triangle2CFmesh::setup() => end\n");
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::unsetup()
{
  LogToScreen (VERBOSE, "Triangle2CFmesh::unsetup()\n");

  m_param2prim.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::configure(OptionMap& cmap, const std::string& prefix)
{
  Converter::configure(cmap, prefix);

  // assign strings read on input.case file to variable transformer object
  m_param2prim.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;

 
  if (ConfigFileReader::isFirstConfig()) {  
   m_param2prim.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                             getProvider(m_param2prim.name())
                             ->create(m_param2prim.name()));
  }

  // configure variable transformer object
  configureDeps(cmap, m_param2prim.ptr().get());
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::convert()
{
  LogToScreen (VERBOSE, "Triangle2CFmesh::convert()\n");

  setMeshData();
  setPhysicsData();
  
  setAddress();

  /// read triangle format file
  LogToScreen(DEBUG_MIN, "Triangle2CFmesh::reading Triangle format\n");
  readTriangleFmt();

  /// make the tansformation from Roe parameter vector variables
  /// to primitive dimensional variables for CFmesh format
  m_param2prim.ptr()->transform();

  /// write CFmesh format
  LogToScreen(DEBUG_MIN, "Triangle2CFmesh::writing CFmesh format\n");
  writeCFmeshFmt();
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::readTriangleFmt()
{
  // dummy variables
  int dummy, iattr, ielem, ivert;
  unsigned ibfac, IFACE, I;
  int in1, in2, n1, n2, IBC;

  // number of boundary faces
  unsigned nbBoundaryfaces;

  /// number of egdes
  unsigned nedge;

  /// reading file
  ifstream file;

  // create Jcycl object
  Jcycl J;

  // set ICLR vector to 0 value
  ICLR.resize(20,0);

  string dummyfile;

  dummyfile = fname->at(0)+".1.node";
  file.open(dummyfile.c_str()); // .node file
  // read number of points
  file >> (*npoin) >> dummy >> dummy >> dummy;
  cout << "ReadTriangleFmt::Opening and reading " << (*npoin);
  cout << " meshpoints from " << fname->at(0) << ".1.node" << endl;
  // resize zroe vector of MeshData pattern
  zroe->resize((*ndofmax) * ((*npoin) + 2 * (*nshmax) * (*npshmax)));
  // assign start pointers for the c_Zroe and c_XY array with the new npoin
  // values
  setAddress();

  // read mesh points status
  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   file >> dummy >> (*c_XY)(0,IPOIN) >> (*c_XY)(1,IPOIN);
   for(unsigned K=0; K<(*ndof); K++) { file >> (*c_Zroe)(K,IPOIN); }
   file >> dummy;
  }
  file.close();
  cout << "ReadTriangleFmt::file node => Done!\n";

  dummyfile = fname->at(0)+".1.ele";
  file.open(dummyfile.c_str()); // .ele file
  // read number of elements
  file >> (*nelem) >> dummy >> dummy;
  cout << "ReadingTriangleFmt::Opening and reading " << (*nelem);
  cout << " triangles from " << fname->at(0) << ".1.ele" << endl;
  // resize arrays whit the new nelem value reas on ele file
  celcel->resize(3,*nelem);
  celnod->resize(3,*nelem);
  // read celnod array
  for(unsigned IELEM=0; IELEM<(*nelem); IELEM++) {
   file >> dummy;
   for(unsigned J=0; J<3; J++)  { file >> (*celnod)(J,IELEM); }
  }
  file.close();
  cout << "ReadTriangleFmt::file ele => Done!\n";

  dummyfile = fname->at(0)+".1.neigh";
  file.open(dummyfile.c_str()); // .neigh file
  // read celcel array
  file >> (*nelem) >> dummy;
  cout << "ReadTriangleFmt::Opening and reading " << (*nelem);
  cout << " neighbours from " << fname->at(0) << ".1.neigh" << endl;
  for(unsigned IELEM=0; IELEM<(*nelem); IELEM++) {
   file >> dummy;
   for(unsigned J=0; J<3; J++)  { file >> (*celcel)(J,IELEM); }
  }
  file.close();
  cout << "ReadTriangleFmt::file neigh => Done!\n";

  dummyfile = fname->at(0)+".1.poly";
  file.open(dummyfile.c_str()); // .poly file
  // read number of faces
  file >> dummy >> dummy >> dummy >> dummy;
  file >> (*nbfac);
  // resize array with the new nbfac value read on poly file
  bndfac->resize(3,(*nbfac) + 2 * (*nshmax) * (*neshmax));
  file.close();
  cout << "ReadTriangleFmt::file poly => Done!\n";

  // fill bndfac values
  // only for element bndfac(0,*) and bndfac(1,*)
  // bndfac(2,*) will be filled after
  ibfac=0;
  for(unsigned IELEM=0; IELEM<(*nelem); IELEM++) {
   for(unsigned I=0; I<3; I++) {
    if((*celcel)(I,IELEM)==-1) {
     (*celcel)(I,IELEM)=0;
     in1 = (*celnod)(J.callJcycl(I+2)-1,IELEM); // c++ indeces start from 0
     in2 = (*celnod)(J.callJcycl(I+3)-1,IELEM); // c++ indeces start from 0
     (*bndfac)(0,ibfac) = IELEM+1; // c++ indeces start from 0
     (*bndfac)(1,ibfac) = I+1; // c++ indeces start from 0
     (*bndfac)(2,ibfac) = -1; // still undefined
     ++ibfac;
    }
   }
  }
  cout << "ReadTriangleFmt::There appear to be " << ibfac << " boundary faces\n";

  nbBoundaryfaces = ibfac;

  dummyfile = fname->at(0)+".1.edge";
  file.open(dummyfile.c_str()); // .edge file
  // read number of edges
  file >> nedge >> iattr;
  cout << "ReadTriangleFmt::Opening and reading " << nedge << " faces from ";
  cout << fname->at(0) << ".1.edge" << endl;
  if (iattr!=1) {
   cout << "ReadTriangleFmt::There should be 1 bndry marker in ";
   cout << fname->at(0) << ".1.edge while there appear to be " << iattr;
   cout << "\n Run Triangle with -e\n";
   exit(1);
  }

  IFACE=0; I=0; ibfac = 0;
  
  while(IFACE<nedge) {
     file >> dummy >> n1 >> n2 >> IBC;
     // if a boundary faces, look for the parent element in bndfac
     if(IBC==0) { ++ICLR.at(IBC); }
     else {
      ++ibfac;
      I=0;
      while(I<(nbBoundaryfaces)) {
       
       ielem = (*bndfac)(0,I);
       ivert = (*bndfac)(1,I);
       
       in1 = (*celnod)(J.callJcycl(ivert+1)-1,ielem-1);//c++ indeces start from 0
       in2 = (*celnod)(J.callJcycl(ivert+2)-1,ielem-1);//c++ indeces start from 0
       
       if ((in1==n1 && in2==n2) || (in1==n2 && in2==n1)) {
        if((*bndfac)(2,I)!=-1) { cout << "ReadTriangleFmt::error => " << I;
                                 cout << " " << (*bndfac)(2,I) << "\n";    
                                 exit(1);                                  }
        (*bndfac)(2,I) = IBC;
        ICLR.at(IBC) = ICLR.at(IBC) + 1;
        goto nine;
       } // if (in1==n1 && in2==n2) || (in1==n2 && in2==n1)
      I++;
      } // while I<(nbBoundaryfaces)
      cout << "ReadTriangleFmt::error => Cannot match boundary face ";
      cout << IFACE << " " << n1 << " " << n2 << " " << IBC << "\n";
      exit(IFACE);
     }
nine:
   IFACE++;
   }
   cout << "ReadTriangle::Boundary faces are " << ibfac << " ";
   cout << nbBoundaryfaces << endl;

   file.close();
   cout << "ReadTriangleFmt::file edge => Done!\n";

  for(unsigned K=0; K<20; K++) {
   if(ICLR.at(K)!=0) {
    cout << "ReadTriangleFmt::Boundary edges colored " << K << " are ";
    cout << ICLR.at(K) << "\n"; }
   if (ibfac!=nbBoundaryfaces) {
    cout << "ReadTriangleFmt::error => No matching boundary faces\n";
    exit(1);
   }
  }
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::writeCFmeshFmt()
{
  unsigned BNDS=0; unsigned BND=0; unsigned IND2=2;
  // @param LIST_STATE = 0 there is not a list of states
  // @param LIST_STATE = 1 there is a list of states
  unsigned  LIST_STATE=1;
  int ip;
  vector <int> np(2);

  // writing file
  ofstream file;

  // create Jcycl object
  Jcycl J;

  // find max value in bndfac(2,*) vector
  int maxNCl = (*bndfac)(2,0);
  for(unsigned IBFAC=0; IBFAC<(*nbfac); IBFAC++) {
   if((*bndfac)(2,IBFAC)>maxNCl) { maxNCl = (*bndfac)(2,IBFAC); }
  }

  for(int IBC=0; IBC<maxNCl; IBC++) {
   if(ICLR.at(IBC+1)>0) { ++BNDS; }
  }

  file.open("cfin.CFmesh");
  file.setf(ios::adjustfield);
  file.precision(20);

  file << "!NB_DIM " << (*ndim) << "\n";
  file << "!NB_EQ " << (*ndof) << "\n";
  file << "!NB_NODES " << (*npoin) << " 0\n";
  file << "!NB_STATES "<< (*npoin) << " 0\n";
  file << "!NB_ELEM " << (*nelem) << "\n";
  file << "!NB_ELEM_TYPES " << "1\n";
  file << "!GEOM_POLYORDER " << "1\n";
  file << "!ELEM_TYPES " << "Triag\n";
  file << "!NB_ELEM_PER_TYPE " << (*nelem) << "\n";
  file << "!NB_NODES_PER_TYPE " << "3\n";
  file << "!NB_STATES_PER_TYPE " << "3\n";
  file << "!LIST_ELEM " << "\n";
  for(unsigned IELEM=0; IELEM<(*nelem); IELEM++) {
   file << (*celnod)(0,IELEM)-1 << " " << (*celnod)(1,IELEM)-1 << " ";
   file << (*celnod)(2,IELEM)-1 << " " << (*celnod)(0,IELEM)-1 << " ";
   file << (*celnod)(1,IELEM)-1 << " " << (*celnod)(2,IELEM)-1 << "\n";
  }

  file << "!NB_TRSs " << BNDS << "\n";

  for(int IBC=0; IBC<maxNCl; IBC++) {
   if(ICLR.at(IBC+1)>0) { ++BND; }
   if ((IBC+1)==10) { file << "!TRS_NAME " << "10\n"; }
   else         { file << "!TRS_NAME " << BND << "\n"; }

   file << "!NB_TRs " << "1\n";
   file << "!NB_GEOM_ENTS " << ICLR.at(IBC+1) << "\n";
   file << "!GEOM_TYPE " << "Face\n";
   file << "!LIST_GEOM_ENT" << "\n";

   for(unsigned j=0; j<(*nbfac); j++) {
    if((*bndfac)(2,j)==(IBC+1)) {
     int ielem = (*bndfac)(0,j);
     int ivert = (*bndfac)(1,j);
     for(unsigned k=0; k<2; k++) {
      ip = (*celnod)(J.callJcycl(ivert+k+1)-1,ielem-1); // c++ indeces start from 0
      np.at(k) = ip-1;
     }
     file << IND2 << " " << IND2 << " ";
     file << np.at(0) << " " << np.at(1) << " ";
     file << np.at(0) << " " << np.at(1) << "\n";
    }
   }
  }

  file << "!LIST_NODE " << "\n";
  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   file << (*c_XY)(0,IPOIN) << " " << (*c_XY)(1,IPOIN) << "\n";
  }


  file << "!LIST_STATE " << LIST_STATE << "\n";
  if(LIST_STATE==1) {
   for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     file << (*c_Zroe)(K,IPOIN) << " ";}
    file << "\n";
   }
  }

  file << "!END " << "\n";
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setAddress()
{
  unsigned start;
  start = 0;
  c_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
  c_Zroe = new Array2D <double> ((*ndof),(*npoin),&zroe->at(start));
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nelem = MeshData::getInstance().getData <unsigned> ("NELEM");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  zroe = MeshData::getInstance().getData< std::vector<double> >("ZROE");
  coor = MeshData::getInstance().getData< std::vector<double> >("COOR");
  celnod = MeshData::getInstance().getData< Array2D<int> >("CELNOD");
  celcel = MeshData::getInstance().getData< Array2D<int> >("CELCEL");
  bndfac = MeshData::getInstance().getData< Array2D<int> >("BNDFAC");
  fname = MeshData::getInstance().getData<vector<string> >("FNAME");
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setPhysicsData()
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
