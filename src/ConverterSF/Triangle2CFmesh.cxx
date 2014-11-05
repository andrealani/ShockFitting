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
  delete XY; delete zroe; delete celnod; delete celcel; delete bndfac;
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
  LogToScreen (INFO, "Triangle2CFmesh::convert()\n");

  setMeshData();
  setPhysicsData();
  
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
  file >> npoin->at(1) >> dummy >> dummy >> dummy;

  // resize zroe vectors of MeshData pattern
  totsize = npoin->at(0) + npoin->at(1) + 4 * (*nshmax) * (*npshmax);
  zroeVect->resize((*ndofmax) * totsize);
  coorVect->resize((*ndim) * totsize);

  // assign start pointers for the zroe and XY arrays
  start = (*ndim) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  XY = new Array2D <double> ((*ndim),
                                  (npoin->at(1) + 2 * (*nshmax) * (*npshmax)),
                                  &coorVect->at(start));
  start = (*ndofmax) * (npoin->at(0) + 2 * (*nshmax) * (*npshmax));
  zroe = new Array2D <double> ((*ndofmax),
                                    (npoin->at(1) + 2 * (*nshmax) * (*npshmax)),
                                    &zroeVect->at(start));

  // read mesh points status
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   file >> dummy >> (*XY)(0,IPOIN) >> (*XY)(1,IPOIN);
   for(unsigned K=0; K<(*ndof); K++) { file >> (*zroe)(K,IPOIN); }
   file >> dummy;
  }
  file.close();

  dummyfile = fname->at(0)+".1.ele";
  file.open(dummyfile.c_str()); // .ele file
  // read number of elements
  file >> nelem->at(1) >> (*nvt) >> dummy;

  // resize arrays whit the new nelem value read on ele file
  totsize = nelem->at(0) + nelem->at(1);
  celcelVect->resize((*nvt) * totsize);
  celnodVect->resize((*nvt) * totsize);
  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  celcel = new Array2D<int> ((*nvt), nelem->at(1), &celcelVect->at(start));

  // read celnod array
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   file >> dummy;
   for(unsigned J=0; J<3; J++)  { file >> (*celnod)(J,IELEM); }
  }
  file.close();

  dummyfile = fname->at(0)+".1.neigh";
  file.open(dummyfile.c_str()); // .neigh file
  // read celcel array
  file >> nelem->at(1) >> dummy;
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   file >> dummy;
   for(unsigned J=0; J<3; J++)  { file >> (*celcel)(J,IELEM); }
  }
  file.close();

  dummyfile = fname->at(0)+".1.poly";
  file.open(dummyfile.c_str()); // .poly file
  // read number of faces
  file >> dummy >> dummy >> dummy >> dummy;
  file >> nbfac->at(1);

  // resize array with the new nbfac value read on poly file
  totsize = nbfac->at(0) + nbfac->at(1) + 4 * (*nshmax) * (*neshmax);
  bndfacVect->resize(3 * totsize);
  start = 3* (nbfac->at(0) + 2 * (*nshmax) * (*neshmax));
  bndfac = new Array2D<int> (3,(nbfac->at(1) + 2 * (*nshmax) * (*neshmax)),
                             &bndfacVect->at(start));

  file.close();

  // fill bndfac values
  // only for element bndfac(0,*) and bndfac(1,*)
  // bndfac(2,*) will be filled after
  ibfac=0;
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
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

  nbBoundaryfaces = ibfac;

  dummyfile = fname->at(0)+".1.edge";
  file.open(dummyfile.c_str()); // .edge file
  // read number of edges
  file >> nedge >> iattr;
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

  file.close();

  if (ibfac!=nbBoundaryfaces) {
   cout << "ReadTriangleFmt::error => No matching boundary faces\n";
   exit(1);
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
  for(unsigned IBFAC=0; IBFAC<nbfac->at(1); IBFAC++) {
   if((*bndfac)(2,IBFAC)>maxNCl) { maxNCl = (*bndfac)(2,IBFAC); }
  }

  for(int IBC=0; IBC<maxNCl; IBC++) {
   if(ICLR.at(IBC+1)>0) { ++BNDS; }
  }

  file.open("cfin.CFmesh");

  file << "!NB_DIM " << setw(1) << (*ndim) << "\n";
  file << "!NB_EQ " << setw(1) << (*ndof) << "\n";
  file << "!NB_NODES " << setw(5) << npoin->at(1) << " 0\n";
  file << "!NB_STATES "<< setw(5) << npoin->at(1) << " 0\n";
  file << "!NB_ELEM " << setw(5) << nelem->at(1) << "\n";
  file << "!NB_ELEM_TYPES 1\n";
  file << "!GEOM_POLYORDER 1\n";
  file << "!SOL_POLYORDER 1\n";
  file << "!ELEM_TYPES Triag\n";
  file << "!NB_ELEM_PER_TYPE " << setw(5) << nelem->at(1) << "\n";
  file << "!NB_NODES_PER_TYPE 3\n";
  file << "!NB_STATES_PER_TYPE 3\n";
  file << "!LIST_ELEM" << "\n";
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   file.setf(ios::right,ios::adjustfield);
   file << " " << setw(10) << (*celnod)(0,IELEM)-1;
   file << " " << setw(10) << (*celnod)(1,IELEM)-1;
   file << " " << setw(10) << (*celnod)(2,IELEM)-1;
   file << " " << setw(10) << (*celnod)(0,IELEM)-1;
   file << " " << setw(10) << (*celnod)(1,IELEM)-1;
   file << " " << setw(10) << (*celnod)(2,IELEM)-1 << "\n";
  }

  file << "!NB_TRSs " << setw(3) << BNDS << "\n";

  for(int IBC=0; IBC<maxNCl; IBC++) {

   if(ICLR.at(IBC+1)>0) { 
    ++BND; 
 
   if ((IBC+1)==10) { file << "!TRS_NAME  10\n"; }

    else            { file << "!TRS_NAME " << setw(3) << BND << "\n"; }

    file << "!NB_TRs 1\n";
    file << "!NB_GEOM_ENTS" << setw(5) << ICLR.at(IBC+1) << "\n";
    file << "!GEOM_TYPE Face\n";
    file << "!LIST_GEOM_ENT" << "\n";

    for(unsigned j=0; j<nbfac->at(1); j++) {
     if((*bndfac)(2,j)==(IBC+1)) {
      int ielem = (*bndfac)(0,j);
      int ivert = (*bndfac)(1,j);
      for(unsigned k=0; k<2; k++) {
       ip = (*celnod)(J.callJcycl(ivert+k+1)-1,ielem-1); // c++ indeces start from 0
       np.at(k) = ip-1;
      }
      file << setw(1) << IND2 << " " << setw(1) << IND2;
      file << " " << setw(10) << np.at(0);
      file << " " << setw(10) << np.at(1);
      file << " " << setw(10) << np.at(0);
      file << " " << setw(10) << np.at(1) << "\n";
     }
    }
   }
  }
  file << "!LIST_NODE" << "\n";
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   file.precision(16);
   file.setf(ios::scientific);
   file << " " << setw(32) << (*XY)(0,IPOIN);
   file << " " << setw(32) << (*XY)(1,IPOIN) << "\n";
  }


  file << "!LIST_STATE " << setw(1) << LIST_STATE << "\n";
  if(LIST_STATE==1) {
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     file.precision(16);
     file.setf(ios::scientific);
     file << setw(32) << (*zroe)(K,IPOIN);}
    file << "\n";
   }
  }

  file << "!END" << "\n";
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  celcelVect = MeshData::getInstance().getData <vector<int> >("CELCEL");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  fname = MeshData::getInstance().getData <vector<string> >("FNAME");
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
