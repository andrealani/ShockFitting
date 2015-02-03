// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/Triangle2CFmeshFreez.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
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
ObjectProvider<Triangle2CFmeshFreez, Converter>
Triangle2CFmeshFreezProv("Triangle2CFmeshFreez");

//--------------------------------------------------------------------------//

Triangle2CFmeshFreez::Triangle2CFmeshFreez(const std::string& objectName) :
  Converter(objectName)
{
  m_boundary = "single";
  addOption("ShockBoundary", &m_boundary,
            "Additional info on the shock boundary: single or splitted");

  m_param2prim.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

Triangle2CFmeshFreez::~Triangle2CFmeshFreez()
{
  delete XY; delete zroe; delete celnod; delete celcel; delete bndfac;
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::setup()
{
  LogToScreen (VERBOSE, "Triangle2CFmeshFreez::setup() => start\n");

  m_param2prim.ptr()->setup();

  LogToScreen (VERBOSE, "Triangle2CFmeshFreez::setup() => end\n");
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::unsetup()
{
  LogToScreen (VERBOSE, "Triangle2CFmeshFreez::unsetup()\n");

  m_param2prim.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::configure(OptionMap& cmap, const std::string& prefix)
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

void Triangle2CFmeshFreez::convert()
{
  LogToScreen (INFO, "Triangle2CFmeshFreez::convert()\n");

  setMeshData();
  setPhysicsData();
  
  if (MeshData::getInstance().getVersion()==("original")) 
  {
   // read triangle format file
   LogToScreen(DEBUG_MIN, "Triangle2CFmeshFreez::reading Triangle format\n");
   readTriangleFmt();
  }

  // make the tansformation from Roe parameter vector variables
  // to primitive dimensional variables for CFmesh format
  m_param2prim.ptr()->transform();

  // write CFmesh format
  LogToScreen(DEBUG_MIN, "Triangle2CFmeshFreez::writing CFmesh format\n");
  writeCFmeshFmt();
  
  // de-allocate dynamic arrays
  freeArray();
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::readTriangleFmt()
{
  // dummy variables
  unsigned dummy, iattr;
  unsigned ibfac, IFACE, I;
  int in1, in2, n1, n2, help, IBC, ielem, ivert;

  // number of boundary faces
  unsigned nbBoundaryfaces;

  /// number of egdes
  unsigned nedge;

  /// reading file
  ifstream file;

  // create Jcycl object
  Jcycl J;

  // set ICLR vector to 0 value
  ICLR->assign(20,0);

  string dummyfile;

  dummyfile = fname->str()+".1.node";

  file.open(dummyfile.c_str()); // .node file
  // read number of points
  file >> npoin->at(1) >> dummy >> dummy >> dummy;

  // resize zroe vectors of MeshData pattern
  totsize = npoin->at(0) + npoin->at(1) + 
            4 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();
  zroeVect->resize(PhysicsInfo::getnbDofMax() * totsize);
  coorVect->resize(PhysicsInfo::getnbDim() * totsize);

  // assign start pointers for the zroe and XY arrays
  start = PhysicsInfo::getnbDim() * 
          (npoin->at(0) + 2 * 
           PhysicsInfo::getnbShMax() *
           PhysicsInfo::getnbShPointsMax());
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 * 
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                             &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() *
          (npoin->at(0) + 2 *
          PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),
                              (npoin->at(1) + 2 *
                               PhysicsInfo::getnbShMax() *
                               PhysicsInfo::getnbShPointsMax()),
                               &zroeVect->at(start));

  // read mesh points state
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   file >> dummy >> (*XY)(0,IPOIN) >> (*XY)(1,IPOIN);
   for(unsigned K=0; K<(*ndof); K++) { file >> (*zroe)(K,IPOIN); }
   file >> dummy;
  }
  file.close();

  dummyfile = fname->str()+".1.ele";

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

  dummyfile = fname->str()+".1.neigh";

  file.open(dummyfile.c_str()); // .neigh file
  // read celcel array
  file >> nelem->at(1) >> dummy;
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   file >> dummy;
   for(unsigned J=0; J<3; J++)  { file >> (*celcel)(J,IELEM); }
  }
  file.close();

  dummyfile = fname->str()+".1.poly";

  file.open(dummyfile.c_str()); // .poly file
  // read number of faces
  file >> dummy >> dummy >> dummy >> dummy;
  file >> nbfac->at(1);

  // resize array with the new nbfac value read on poly file
  totsize = nbfac->at(0) + nbfac->at(1) + 
            4 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax();
  bndfacVect->resize(3 * totsize);
  start = 3 * (nbfac->at(0) + 
          2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax());
  bndfac = new Array2D<int> (3,(nbfac->at(1) +
                             2 * PhysicsInfo::getnbShMax() *
                             PhysicsInfo::getnbShEdgesMax()),
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

  dummyfile = fname->str()+".1.edge";

  file.open(dummyfile.c_str()); // .edge file
  // read number of edges
  file >> nedge >> iattr;
  if (iattr!=1) {
   cout << "ReadTriangleFmt::There should be 1 bndry marker in ";
   cout << fname->str() << ".1.edge while there appear to be " << iattr;
   cout << "\n Run Triangle with -e\n";
   exit(1);
  }

  IFACE=0; I=0; ibfac = 0;
  
  while(IFACE<nedge) {
     file >> dummy >> n1 >> n2 >> IBC;
     // if a boundary faces, look for the parent element in bndfac
     if(IBC==0 || IBC==999) { help = ICLR->at(0);
                              ICLR->at(0) = help +1; }
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
        help = ICLR->at(IBC);
        ICLR->at(IBC) = help + 1;
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

void Triangle2CFmeshFreez::writeCFmeshFmt()
{
  unsigned BNDS=0; unsigned BND=0; unsigned IND2=2;
  // @param LIST_STATE = 0 there is not a list of states
  // @param LIST_STATE = 1 there is a list of states
  unsigned  LIST_STATE=1;
  int ip;
  vector <unsigned> np(2);

  // writing file
  FILE* cfin;

  // create Jcycl object
  Jcycl J;

  // allocate the arrays if the optimized version
  if(MeshData::getInstance().getVersion()=="optimized")
  {
   // assign start pointers for the zroe and XY arrays
   start = PhysicsInfo::getnbDim() *
           (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() *
           PhysicsInfo::getnbShPointsMax());
   XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                             &coorVect->at(start));
   start = PhysicsInfo::getnbDofMax() *
           (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
   zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),
                               (npoin->at(1) + 2 *
                                PhysicsInfo::getnbShMax() *
                                PhysicsInfo::getnbShPointsMax()),
                                &zroeVect->at(start));
  // assign starting pointers for the celcel and celnod arrays
   start = (*nvt) * nelem->at(0);
   celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
   celcel = new Array2D<int> ((*nvt), nelem->at(1), &celcelVect->at(start));
  // assign the starting pointers for the bndfac array
   start = 3 * (nbfac->at(0) +
          2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax());
   bndfac = new Array2D<int> (3,(nbfac->at(1) +
                              2 * PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShEdgesMax()),
                              &bndfacVect->at(start));
  }
 
  // find max value in bndfac(2,*) vector
  int maxNCl = (*bndfac)(2,0);
  for(unsigned IBFAC=0; IBFAC<nbfac->at(1); IBFAC++) {
   if((*bndfac)(2,IBFAC)>maxNCl) { maxNCl = (*bndfac)(2,IBFAC); }
  }

  for(int IBC=0; IBC<maxNCl; IBC++) {
   if(ICLR->at(IBC+1)>0) { ++BNDS; }
  }

  // compute number of shock elements
  int nbSh = ICLR->at(10) + 2;  

  unsigned minSh=npoin->at(1); unsigned maxSh=0;
  for(unsigned IFACE=0; IFACE<nbfac->at(0); IFACE++) {
   if((*bndfac)(2,IFACE)==10) {
    int elem = (*bndfac)(0,IFACE);
    int vert = (*bndfac)(1,IFACE);
    for(unsigned k=0; k<2; k++) {
     ip = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1); // c++ indeces start from 0
     np.at(k) = ip-1;
     if (np.at(k) < minSh) { minSh = np.at(k); }
     if (np.at(k) > maxSh) { maxSh = np.at(k); }
    } // for k<2
   } // if (*bndfac)(2,IFACE)==10
  } // for IFACE<nbfac->at(0)


  cfin = fopen("cfin.CFmesh", "w");

  fprintf(cfin,"%s %1u","!NB_DIM",PhysicsInfo::getnbDim());
  fprintf(cfin,"%s %1u","\n!NB_EQ",(*ndof));
  fprintf(cfin,"%s %5u %s","\n!NB_NODES",npoin->at(1) ,"0\n");
  fprintf(cfin,"%s %5u %s","!NB_STATES",npoin->at(1),"0\n");
  fprintf(cfin,"%s %5u","!NB_ELEM",nelem->at(1));
  fprintf(cfin,"%s","\n!NB_ELEM_TYPES 1\n");
  fprintf(cfin,"%s","!GEOM_POLYORDER 1\n");
  fprintf(cfin,"%s","!SOL_POLYORDER 1\n");
  fprintf(cfin,"%s","!ELEM_TYPES Triag\n");
  fprintf(cfin,"%s %5u","!NB_ELEM_PER_TYPE",nelem->at(1));
  fprintf(cfin,"%s","\n!NB_NODES_PER_TYPE 3\n");
  fprintf(cfin,"%s","!NB_STATES_PER_TYPE 3\n");
  fprintf(cfin,"%s","!LIST_ELEM");
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   fprintf(cfin,"%s %10i","\n",(*celnod)(0,IELEM)-1);
   fprintf(cfin,"%11i",(*celnod)(1,IELEM)-1);
   fprintf(cfin,"%11i",(*celnod)(2,IELEM)-1);
   fprintf(cfin,"%11i",(*celnod)(0,IELEM)-1);
   fprintf(cfin,"%11i",(*celnod)(1,IELEM)-1);
   fprintf(cfin,"%11i",(*celnod)(2,IELEM)-1);
  }

  if (m_boundary == "single") {
   fprintf(cfin,"%s %3u","\n!NB_TRSs",BNDS); 
  }

  else if (m_boundary == "splitted") {
   fprintf(cfin,"%s %3u","\n!NB_TRSs",BNDS+1);
  }

  for(int IBC=0; IBC<maxNCl; IBC++) {

   if(ICLR->at(IBC+1)>0) { 
    ++BND; 
    if((IBC+1)==10) {
     if (m_boundary == "single") {
      fprintf(cfin,"%s","\n!TRS_NAME  10\n");
      fprintf(cfin,"%s","!NB_TRs 1\n");
      fprintf(cfin,"%s %4u","!NB_GEOM_ENTS",ICLR->at(IBC+1));
      fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
      fprintf(cfin,"%s","!LIST_GEOM_ENT");

      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         ip = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1); // c++ indeces start from 0
         np.at(k) = ip-1;
        }
        fprintf(cfin,"%s","\n");
        fprintf(cfin,"%1i %1i",IND2,IND2);
        fprintf(cfin,"%11i %10i %10i %10i",np.at(0),np.at(1),np.at(0),np.at(1));
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)
     } // if m_boundary = single

     if (m_boundary == "splitted") {
      
      // Supersonic boundary
      fprintf(cfin,"%s","\n!TRS_NAME  InnerSup\n");
      fprintf(cfin,"%s","!NB_TRs 1\n");
      fprintf(cfin,"%s %4u","!NB_GEOM_ENTS",nbSh/2-1);
      fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
      fprintf(cfin,"%s","!LIST_GEOM_ENT");

      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         ip = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1); // c++ indeces start from 0
         np.at(k) = ip-1;
        }
        if ((np.at(0) >= minSh) && (np.at(0) <  (minSh+nbSh/2)) &&
            (np.at(1) >= minSh) && (np.at(1) <  (minSh+nbSh/2))) {
        fprintf(cfin,"%s","\n");
        fprintf(cfin,"%1i %1i",IND2,IND2);
        fprintf(cfin,"%11i %10i %10i %10i",np.at(0),np.at(1),np.at(0),np.at(1));
        } // if np conditions
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

      // Subsonic boundary
      fprintf(cfin,"%s","\n!TRS_NAME  InnerSub\n");
      fprintf(cfin,"%s","!NB_TRs 1\n");
      fprintf(cfin,"%s %4u","!NB_GEOM_ENTS",nbSh/2-1);
      fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
      fprintf(cfin,"%s","!LIST_GEOM_ENT");

      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         ip = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1); // c++ indeces start from 0
         np.at(k) = ip-1;
        }
        if ((np.at(0) > (maxSh-nbSh/2)) && (np.at(0) <= maxSh) &&
            (np.at(1) > (maxSh-nbSh/2)) && (np.at(1) <= maxSh)) {
        fprintf(cfin,"%s","\n");
        fprintf(cfin,"%1i %1i",IND2,IND2);
        fprintf(cfin,"%11i %10i %10i %10i",np.at(0),np.at(1),np.at(0),np.at(1));
        } // if np conditions
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

     } // if m_boundary = splitted
    } // if IBC==10

    else            { 
     fprintf(cfin,"%s %3u","\n!TRS_NAME",BND);
     fprintf(cfin,"%s","\n!NB_TRs 1\n");
     fprintf(cfin,"%s %4i","!NB_GEOM_ENTS",ICLR->at(IBC+1));
     fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
     fprintf(cfin,"%s","!LIST_GEOM_ENT");

     for(unsigned j=0; j<nbfac->at(1); j++) {
      if((*bndfac)(2,j)==(IBC+1)) {
       int ielem = (*bndfac)(0,j);
       int ivert = (*bndfac)(1,j);
       for(unsigned k=0; k<2; k++) {
        ip = (*celnod)(J.callJcycl(ivert+k+1)-1,ielem-1); // c++ indeces start from 0
        np.at(k) = ip-1;
       }
       fprintf(cfin,"%s","\n");
       fprintf(cfin,"%1i %1i",IND2,IND2);
       fprintf(cfin,"%11i %10i %10i %10i",np.at(0),np.at(1),np.at(0),np.at(1));
      } // if (*bndfac)(2,j)==(IBC+1)
     } // for j<nbfac->at(1)
    } // else (ICLR(IBC+1)!=10)

   } // if ICLR.at(IBC+1)>0
  } // for IBC<maxNCl


  fprintf(cfin,"%s","\n!LIST_NODE\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
    fprintf(cfin,"%s"," ");
    fprintf(cfin,"%32.16E",(*XY)(IA,IPOIN)); }
   fprintf(cfin,"%s","\n");
  }

  fprintf(cfin,"%s %1u","!LIST_STATE",LIST_STATE);
  fprintf(cfin,"%s","\n");
  if(LIST_STATE==1) {
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     fprintf(cfin,"%s"," ");
     fprintf(cfin,"%32.16E",(*zroe)(K,IPOIN));}
    fprintf(cfin,"%s","\n");
   }
  }

  fprintf(cfin,"%s","!END\n");

  fclose(cfin);
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::freeArray()
{
  delete XY; delete zroe;
  delete bndfac; delete celnod; delete celcel;
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::setMeshData()
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
  fname = MeshData::getInstance().getData <stringstream>("FNAME");
  ICLR = MeshData::getInstance().getData <vector<int> >("ICLR");
}

//----------------------------------------------------------------------------//

void Triangle2CFmeshFreez::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
