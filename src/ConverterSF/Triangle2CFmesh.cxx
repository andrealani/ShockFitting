// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/Triangle2CFmesh.hh"
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
ObjectProvider<Triangle2CFmesh, Converter>
triangle2CFmeshProv("Triangle2CFmesh");

//--------------------------------------------------------------------------//

Triangle2CFmesh::Triangle2CFmesh(const std::string& objectName) :
  Converter(objectName)
{
  m_boundary = "single";
  addOption("ShockBoundary", &m_boundary,
            "Additional info on the shock boundary: single or splitted");

  m_param2prim.name() = "dummyVariableTransformer";
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
  LogToScreen (INFO, "Triangle2CFmesh::convert()\n");

  setMeshData();
  setPhysicsData();

  if (MeshData::getInstance().getVersion()==("original")) 
  {
   // read triangle format file
   LogToScreen(DEBUG_MIN, "Triangle2CFmesh::reading Triangle format\n");
   readTriangleFmt();
  }

  // make the tansformation from Roe parameter vector variables
  // to primitive dimensional variables for CFmesh format
  m_param2prim.ptr()->transform();

  // write CFmesh format
  LogToScreen(DEBUG_MIN, "Triangle2CFmesh::writing CFmesh format\n");
  writeCFmeshFmt();

  // store CFmesh data to exchange data with CF without I/O
  LogToScreen(DEBUG_MIN, "Triangle2CFmesh::storing CFmesh data\n");
//  storeCFmeshData();

  // de-allocate dynamic arrays
  freeArray();
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::readTriangleFmt()
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
     if(IBC==0) { help = ICLR->at(IBC);
                  ICLR->at(IBC) = help +1; }
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

void Triangle2CFmesh::writeCFmeshFmt()
{
  unsigned BND=0; unsigned IND2=2;
  // @param LIST_STATE = 0 there is not a list of states
  // @param LIST_STATE = 1 there is a list of states
  unsigned  LIST_STATE = 1;
  int ip;
  vector <unsigned> np(2);

  // variables used to write on the farfield tecplot file
  vector<string> rhostr((*nsp));

  vector<int>elemVector;
  vector<int>elemVector_noduplicate;
  vector<int>elemVector_unitIndex;
  vector<int>elemVector_sort;
  Array2D <int> array_elemVector_sort;

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

  BNDS = 0;
  for(int IBC=0; IBC<maxNCl; IBC++) {
   if(ICLR->at(IBC+1)>0) { ++BNDS; }
  }

  // compute number of shock elements
  int nbSh = ICLR->at(10) + 2;  

  unsigned minSh=npoin->at(1); unsigned maxSh=0;
//  for(unsigned IFACE=0; IFACE<nbfac->at(0); IFACE++) {
  for(unsigned IFACE=0; IFACE<nbfac->at(1); IFACE++) {
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
  fprintf(cfin,"%s","!LIST_ELEM\n");
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
   boundaryNames->resize(BNDS);
  }

  else if (m_boundary == "splitted") {
   fprintf(cfin,"%s %3u","\n!NB_TRSs",BNDS+1);
   boundaryNames->resize(BNDS+1);
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

      // assign the boundary name to the corresponding MeshData vector
      boundaryNames->at(BND-1)="10";

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

      // assign the boundary name to the corresponding MeshData vector
      boundaryNames->at(BND-1) = "InnerSup";

      // Subsonic boundary
      fprintf(cfin,"%s","\n!TRS_NAME  InnerSub\n");
      fprintf(cfin,"%s","!NB_TRs 1\n");
      fprintf(cfin,"%s %4u","!NB_GEOM_ENTS",nbSh/2-1);
      fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
      fprintf(cfin,"%s","!LIST_GEOM_ENT");

      elemVector.resize((nbSh/2-1)*2);
      elemVector_noduplicate.resize((nbSh/2-1)*2);
      elemVector_unitIndex.resize((nbSh/2-1)*2);

      unsigned h=0; unsigned duplicate=0;
      
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
        for(unsigned k=0; k<2; k++) {elemVector.at(h+k)= np.at(k);}
        h=h+2;
        } // if np conditions
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

      // delete the duplicate of the cell nodes IDs from the elemVector
      // by defining an elemVector_noduplicate vector
      // and count the number of not-duplicated nodes
      h=0;
      for(unsigned j=0;j<elemVector.size();j++){
       for(unsigned k=j+1;k<elemVector.size();k++) {
        if(elemVector.at(j)==elemVector.at(k)) { duplicate++; }
       }
       if(duplicate==0) { elemVector_noduplicate.at(h)=elemVector.at(j);
                          h++;
                         }
       else { duplicate=0; }
      }

      elemVector_noduplicate.resize(h);
      elemVector_sort.resize(h);
      array_elemVector_sort.resize(2,h);

      // define a new vector sorting the elements inside elemVector_noduplicate
      elemVector_sort=elemVector_noduplicate;

      // sort the elemVector_sort vector
      sort(elemVector_sort.begin(),elemVector_sort.end());

      // since the tecplot cell nodes for each TRS are numbered 
      // starting from 1:nbElem(TRS), the array_elemVector is defined to 
      //  link the SF cell node IDs to the new IDs defined inside the tecplot file
      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       array_elemVector_sort(0,j)=elemVector_sort.at(j);
       array_elemVector_sort(1,j)=j+1;
      }

      // assign the new IDs to the elemVector_uniIndex vector
      for(unsigned j=0;j<elemVector.size(); j++) {
       for(unsigned k=0;k<array_elemVector_sort.getnCols();k++) {
        if(elemVector.at(j)==array_elemVector_sort(0,k)) {
         elemVector_unitIndex.at(j)=array_elemVector_sort(1,k); }
       }
      }

      // variable storing curvilinear coordinates of the
      // innersub grid points
      double r;

      // output file storing info on the subsonic boundary patch
      FILE* fileFarFieldBC = fopen("FarFieldBc.dat", "w");
      FILE* printFarFieldBC = fopen("FarFieldBc.plt", "w");

      fprintf(fileFarFieldBC,"%s","VARIABLES  =  \"y\" ");
      fprintf(printFarFieldBC,"%s","TITLE = Boundary data\n");
      fprintf(printFarFieldBC,"%s","VARIABLES  =  \"x0\" \"x1\" ");

      if(m_modelTransf=="TCneq") {
      for(unsigned isp=0; isp<(*nsp);isp++) {
       fprintf(fileFarFieldBC,"%s",rhostr.at(isp).c_str());
      }
      fprintf(fileFarFieldBC,"%s","\"u\" \"v\" \"T\" \"Tv0\" \n");
      fprintf(printFarFieldBC,"%s","\"u\" \"v\" \"T\" \"Tv0\" \n");
      }
      else if (m_modelTransf=="Pg") {
       fprintf(fileFarFieldBC,"%s","\"p\" \"u\" \"v\" \"T\"\n");
       fprintf(printFarFieldBC,"%s","\"p\" \"u\" \"v\" \"T\"\n");
      }
      else { cout << m_modelTransf << " model not implemented\n"; exit(1); }

      fprintf(fileFarFieldBC,"%u %s",elemVector_noduplicate.size(),"\n");

      fprintf(printFarFieldBC,"%s %u","ZONE N=",elemVector_noduplicate.size());
      fprintf(printFarFieldBC,"%s %u %s",", T=\"InnerSub, TR 0\", E=",nbSh/2-1,"," );
      fprintf(printFarFieldBC,"%s","F=FEPOINT, ET=LINESEG, SOLUTIONTIME=0\n");
 
      for(unsigned j=0;j<elemVector_noduplicate.size();j++) {
       r = 0;
       for(unsigned k=0;k<PhysicsInfo::getnbDim();k++) {
        r = r + pow((*XY)(k,elemVector_noduplicate.at(j)),2);
        fprintf(printFarFieldBC,"%20.16E %s",(*XY)(k,elemVector_sort.at(j))," ");
       }
         fprintf(fileFarFieldBC,"%28.16e",(*XY)(1,elemVector_noduplicate.at(j)));
       for(unsigned k=0;k<(*ndof);k++) {
        fprintf(fileFarFieldBC,"%30.16e",(*zroe)(k,elemVector_noduplicate.at(j)));
        fprintf(printFarFieldBC,"%30.16e",(*zroe)(k,elemVector_sort.at(j)));
       }
       fprintf(fileFarFieldBC,"%s","\n");
       fprintf(printFarFieldBC,"%s","\n");
      }

      for(unsigned j=0;j<elemVector.size()-1;j=j+2) {
       fprintf(printFarFieldBC,"%11i",elemVector_unitIndex.at(j));
       fprintf(printFarFieldBC,"%11i %s",elemVector_unitIndex.at(j+1),"\n");
      }


      // close file storing InnerSub info
      fclose(fileFarFieldBC);
      fclose(printFarFieldBC);

      // assign the boundary name to the corresponding MeshData vector
      boundaryNames->at(BND) = "InnerSub";

     } // if m_boundary = splitted
    } // if IBC==10

    else            { 
     fprintf(cfin,"%s %3u","\n!TRS_NAME",BND);
     fprintf(cfin,"%s","\n!NB_TRs 1\n");
     fprintf(cfin,"%s %4i","!NB_GEOM_ENTS",ICLR->at(IBC+1));
     fprintf(cfin,"%s","\n!GEOM_TYPE Face\n");
     fprintf(cfin,"%s","!LIST_GEOM_ENT");

     // assign the boundary name to the corresponding MeshData vector
     stringstream dummyBnd;
     dummyBnd << BND;
     boundaryNames->at(BND-1) = dummyBnd.str();

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
    fprintf(cfin,"%32.16E",(*XY)(IA,IPOIN));
   }
   fprintf(cfin,"%s","\n");
  }

  fprintf(cfin,"%s %1u","!LIST_STATE",LIST_STATE);
  fprintf(cfin,"%s","\n");
  if(LIST_STATE==1) {
   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    for(unsigned K=0; K<(*ndof); K++) {
     fprintf(cfin,"%s"," ");
     fprintf(cfin,"%32.16E",(*zroe)(K,IPOIN));
    }
    fprintf(cfin,"%s","\n");
   }
  }

  fprintf(cfin,"%s","!END\n");

  fclose(cfin);
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::storeCFmeshData()
{
  unsigned BND=0;
  int ip;
  vector <unsigned> np(2);

  // dummy ostringstream variable used to store the boundary patches names
  ostringstream iBndName;

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

  BNDS = 0;
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

  // resize array of Connectivity objects
  resizeConnectivityArray();

  unsigned IBNODE=0;
  unsigned m_IBNODE = 0;
  for(int IBC=0; IBC<maxNCl; IBC++) {

  if(ICLR->at(IBC+1)>0) { 
   ++BND;

   if((IBC+1)==10) {
    if (m_boundary == "single") {

     for(unsigned j=0; j<nbfac->at(1); j++) {
      if((*bndfac)(2,j)==(IBC+1)) {
       int elem = (*bndfac)(0,j);
       int vert = (*bndfac)(1,j);
       for(unsigned k=0; k<2; k++) {
        ip = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1); // c++ indeces start from 0
        np.at(k) = ip-1;
       }
       // store the array of the BoundaryConnectivity object 

       // !! it should be more generic and not with only two nodes for each boundary face
       for(unsigned k=0; k<2; k++) {
        // face-node connectivity
        boundaryNodes->at(IBNODE+k)=np.at(k);
        // pointer to the start of the corresponding element
        boundaryPtr->at(IBNODE+k)=IBNODE+k;
       }

       IBNODE=IBNODE+2;

      } // if (*bndfac)(2,j)==(IBC+1)
     } // for j<nbfac->at(1)

     // store the array of the BoundaryConnectivity object
     // array of size @see BNDS x 2, storing consecutively
     // (1) the number of elements (faces) in each boundary
     boundaryInfo->at((BND-1)*2)=ICLR->at(IBC+1);

     // (2) a pointer to the first element in each boundary
     // !! it should be more generic as IBC*nbFace (not only 2)
     boundaryInfo->at((BND-1)*2+1)=boundaryPtr->at(m_IBNODE);
     
     // store the name of the boundary patch
     iBndName << 10;
     boundaryNames->at((BND-1)) = iBndName.str();
     iBndName.str(string());

    } // if m_boundary = single

    if (m_boundary == "splitted") {
     
     // Supersonic boundary
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

       // store the array of the BoundaryConnectivity object

       // !! it should be more generic and not with only two nodes for each boundary face
       for(unsigned k=0; k<2; k++) {
        // face-node connectivity 
        boundaryNodes->at(IBNODE+k)=np.at(k);
        // pointer to the start of the corresponding element 
        boundaryPtr->at(IBNODE+k)=IBNODE+k;
       }

       IBNODE=IBNODE+2;

       } // if np conditions
      } // if (*bndfac)(2,j)==(IBC+1)
     } // for j<nbfac->at(1)

     // store the array of the BoundaryConnectivity object
     // array of size @see BNDS x 2, storing consecutively
     // (1) the number of elements (faces) in each boundary
     boundaryInfo->at((BND-1)*2)=nbSh/2-1;
     // (2) a pointer to the first element in each boundary
     // !! it should be more generic as IBC*nbFace (not only 2)
     boundaryInfo->at((BND-1)*2+1)=boundaryPtr->at(m_IBNODE);

     // store the name of the boundary patch
     iBndName << "InnerSup";
     boundaryNames->at((BND-1)) = iBndName.str();
     iBndName.str(string());

     m_IBNODE = m_IBNODE + boundaryInfo->at((BND-1)*2);

     // Subsonic boundary

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

       // store the array of the BoundaryConnectivity object

       // !! it should be more generic and not with only two nodes for each boundary face
       for(unsigned k=0; k<2; k++) {
        // face-node connectivity 
        boundaryNodes->at(IBNODE+k)=np.at(k);
        // pointer to the start of the corresponding element 
        boundaryPtr->at(IBNODE+k)=IBNODE+k;
       }

       IBNODE=IBNODE+2;

       } // if np conditions
      } // if (*bndfac)(2,j)==(IBC+1)
     } // for j<nbfac->at(1)

     // store the array of the BoundaryConnectivity object
     // array of size @see BNDS x 2, storing consecutively
     // (1) the number of elements (faces) in each boundary
     boundaryInfo->at(BND*2)=nbSh/2-1;
     // (2) a pointer to the first element in each boundary
     // !! it should be more generic as IBC*nbFace (not only 2)
     boundaryInfo->at(BND*2+1)=boundaryPtr->at(m_IBNODE);

     // store the name of the boundary patch
     iBndName << "InnerSub";
     boundaryNames->at(BND) = iBndName.str();
     iBndName.str(string());

    } // if m_boundary = splitted
   } // if IBC==10

   else            { 

    for(unsigned j=0; j<nbfac->at(1); j++) {
     if((*bndfac)(2,j)==(IBC+1)) {
      int ielem = (*bndfac)(0,j);
      int ivert = (*bndfac)(1,j);
      for(unsigned k=0; k<2; k++) {
       ip = (*celnod)(J.callJcycl(ivert+k+1)-1,ielem-1); // c++ indeces start from 0
       np.at(k) = ip-1;
      }

      // store the array of the BoundaryConnectivity object

      // !! it should be more generic and not with only two nodes for each boundary face
      for(unsigned k=0; k<2; k++) {
       // face-node connectivity 
       boundaryNodes->at(IBNODE+k)=np.at(k);
       // pointer to the start of the corresponding element 
       boundaryPtr->at(IBNODE+k)=IBNODE+k;
      }

      IBNODE=IBNODE+2;

     } // if (*bndfac)(2,j)==(IBC+1)
    } // for j<nbfac->at(1)

    // store the array of the BoundaryConnectivity object
    // array of size @see BNDS x 2, storing consecutively
    // (1) the number of elements (faces) in each boundary
    boundaryInfo->at((BND-1)*2)=ICLR->at(IBC+1);

    // (2) a pointer to the first element in each boundary
    // !! it should be more generic as IBC*nbFace (not only 2)
    boundaryInfo->at((BND-1)*2+1)=boundaryPtr->at(m_IBNODE);

    // store the name of the boundary patch
    iBndName << BND;
    boundaryNames->at((BND-1)) = iBndName.str();
    iBndName.str(string());

    } // else (ICLR(IBC+1)!=10)

    m_IBNODE = m_IBNODE + boundaryInfo->at((BND-1)*2);

   } // if ICLR.at(IBC+1)>0
  } // for IBC<maxNCl

  // store element-node connectivity
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   // array storing pointers to the beginning of each element
   elementPtr->at(IELEM)=IELEM*(*nvt)+0;
   // since c++ array start from 0, the element-node connectivity
   // ID are decreased of 1
   // (it does not influence the SF connecitvity data because celnod
   // will be overwritten in CFmesh2Triangle
   for(unsigned IVERT=0;IVERT<(*nvt);IVERT++) {
    (*celnod)(IVERT,IELEM)=(*celnod)(IVERT,IELEM)-1;
   }
  }

  // assign the CF connectivity
  storeConnectivity();
  storeField();
}

//----------------------------------------------------------------------------//

void  Triangle2CFmesh::resizeConnectivityArray()
{
  if(m_boundary=="single") { boundaryInfo->resize(BNDS*2);
                             boundaryNames->resize(BNDS);}
  if(m_boundary=="splitted") { boundaryInfo->resize((BNDS+1)*2);
                               boundaryNames->resize(BNDS+1);}
  // !! it should be more generic and not with only two nodes for each boundary face
  boundaryNodes->resize(nbfac->at(1)*2);
  boundaryPtr->resize(nbfac->at(1)*2);
  elementPtr->resize(nelem->at(1));
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::storeConnectivity()
{
  // @param PhysicsInfo::getnbDim() spatial dimension (2 or 3)
  // @param nbBnds                  number of individual boundary patches
  // @param boundaryInfo            array of size @see BNDS x 2, storing consecutively
  //                                (1) the number of elements (faces) in each boundary
  //                                (2) a pointer to the first element in each boundary
  // @param boundaryNodes           element-node connectivity
  // @param boundaryPtr             pointer to the start of the corresponding element
  inbConnectivity->reset(PhysicsInfo::getnbDim(),
                         boundaryInfo->size()/2,
                         &boundaryInfo->at(0),
                         &boundaryNodes->at(0),
                         &boundaryPtr->at(0),
                         &boundaryNames->at(0));

  // @param nelem->at(1)    number of elements
  // @param celnodVect      element-node connectivity
  // @param elementPtr      pointer to the start of the corresponding element
  inConnectivity->reset(nelem->at(1),
                        &celnodVect->at((*nvt)*nelem->at(0)),
                        &elementPtr->at(0));
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::storeField()
{
  unsigned start;

  // @param npoin           number of degrees of freedom
  // @param inConnectivity  inner connectivity 
  // @param zroeVect        mesh points state
  start = PhysicsInfo::getnbDofMax() *
          (npoin->at(0) + 2 *
          PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  inStateField->reset(npoin->at(1),
                      (*ndof), // stride for dof freedom (e.g.nb of equations) 
                      (*inConnectivity),
                      &zroeVect->at(start));

  // @param npoin             number of degrees of freedom
  // @param inConnectivity    inner connectivity 
  // @param coorVect          mesh points coordinates
  start = PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() *
           PhysicsInfo::getnbShPointsMax());
  inCoordinatesField->reset(npoin->at(1),
                            (*ndof),// stride for dof (e.g.nb of equations) 
                            (*inConnectivity),
                            &coorVect->at(start));
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::freeArray()
{
  delete XY; delete zroe;
  delete celnod; delete celcel; delete bndfac;
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
  fname = MeshData::getInstance().getData <stringstream>("FNAME");
  ICLR = MeshData::getInstance().getData <vector<int> >("ICLR");
  boundaryNodes = MeshData::getInstance().getData <vector<int> >("CF_boundaryNodes");
  boundaryNames =
   MeshData::getInstance().getData <vector<string> >("CF_boundaryNames");
  boundaryInfo = MeshData::getInstance().getData <vector<int> >("CF_boundaryInfo");
  boundaryPtr = MeshData::getInstance().getData <vector<int> >("CF_boundaryPtr");  
  elementPtr = MeshData::getInstance().getData <vector<int> >("CF_elementPtr"); 
  inConnectivity = MeshData::getInstance().getData <Connectivity> ("inConn");
  inbConnectivity = 
    MeshData::getInstance().getData <BoundaryConnectivity> ("inbConn");
  inStateField = MeshData::getInstance().getData <Field> ("inStateField");
  inCoordinatesField = MeshData::getInstance().getData <Field> ("inNodeField");
}

//----------------------------------------------------------------------------//

void Triangle2CFmesh::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
