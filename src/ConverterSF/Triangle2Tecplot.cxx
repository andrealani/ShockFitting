// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/Triangle2Tecplot.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"
#include "MathTools/Isortrx.hh"
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
ObjectProvider<Triangle2Tecplot, Converter>
Triangle2TecplotProv("Triangle2Tecplot");

//--------------------------------------------------------------------------//

Triangle2Tecplot::Triangle2Tecplot(const std::string& objectName) :
  Converter(objectName)
{
  m_boundary = "single";
  addOption("ShockBoundary", &m_boundary,
            "Additional info on the shock boundary: single or splitted");

  m_param2prim.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

Triangle2Tecplot::~Triangle2Tecplot()
{
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::setup()
{
  LogToScreen (VERBOSE, "Triangle2Tecplot::setup() => start\n");

  m_param2prim.ptr()->setup();

  LogToScreen (VERBOSE, "Triangle2Tecplot::setup() => end\n");
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::unsetup()
{
  LogToScreen (VERBOSE, "Triangle2Tecplot::unsetup()\n");

  m_param2prim.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::configure(OptionMap& cmap, const std::string& prefix)
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

void Triangle2Tecplot::convert()
{
  LogToScreen (INFO, "Triangle2Tecplot::convert()\n");

  setMeshData();
  setPhysicsData();
  
  if (MeshData::getInstance().getVersion()==("original")) 
  {
   // read triangle format file
   LogToScreen(DEBUG_MIN, "Triangle2Tecplot::reading Triangle format\n");
   readTriangleFmt();
  }

  // make the tansformation from Roe parameter vector variables
  // to primitive dimensional variables for Tecplot format
  m_param2prim.ptr()->transform();

  // write Tecplot format
  LogToScreen(DEBUG_MIN, "Triangle2Tecplot::writing Tecplot format\n");
  writeTecplotFmt();

  // de-allocate dynamic arrays
  freeArray();
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::readTriangleFmt()
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

  nodcod->resize(npoin->at(1) + 2 * PhysicsInfo::getnbShMax() *
                 PhysicsInfo::getnbShPointsMax());

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
   file >> nodcod->at(IPOIN);
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

void Triangle2Tecplot::writeTecplotFmt()
{
  unsigned BNDS=0; unsigned BND=0; int ip;
  unsigned h, duplicate;

  vector<int>elemVector;
  vector<int>elemVector_noduplicate;
  vector<int>elemVector_unitIndex;
  vector<int>elemVector_sort;
  Array2D <int> array_elemVector_sort;

  // variables used to write on the tecplot files
  vector<string> rhostr((*nsp));
  stringstream nbNodestr;
  stringstream nbElemstr;
  stringstream Tstr;

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

  
  // assign the string of the rho variables
  for(unsigned isp=0;isp<(*nsp);isp++) {
   stringstream ispstr;
   ispstr << isp;
   rhostr.at(isp)="\"rho"+ispstr.str()+"\"";
  }

  // assign the string for the number of nodes
  nbNodestr << "N=" << npoin->at(1);

  // assign the string for the number of elements
  nbElemstr << "E=" << nelem->at(1);

  // write the .plt file
  cfin = fopen("cfin.plt", "w");

  fprintf(cfin,"%s","TITLE      =  Unstructured grid data\n");
  fprintf(cfin,"%s","VARIABLES  =  \"x0\" \"x1\" ");
  for(unsigned isp=0; isp<(*nsp);isp++) {
   fprintf(cfin,"%s %s",rhostr.at(isp).c_str(), " ");
  }
  fprintf(cfin,"%s","\"u\" \"v\" \"T\" \"Tv0\" \n");
  fprintf(cfin,"%s %s","ZONE   T=\"P0 ZONE0 Triag\",",nbNodestr.str().c_str());
  fprintf(cfin,"%s %s %s",",",nbElemstr.str().c_str(),", F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) { 
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
    fprintf(cfin,"%20.16E %s",(*XY)(IA,IPOIN), " ");}
   for(unsigned K=0; K<(*ndof); K++) {
    fprintf(cfin,"%20.16E %s",(*zroe)(K,IPOIN)," ");}  
  fprintf(cfin,"%s","\n");
  } 

  for(unsigned IELEM=0;IELEM<nelem->at(1);IELEM++) {
   for(unsigned K=0;K<3; K++){
    fprintf(cfin,"%11i %s",(*celnod)(K,IELEM)," ");
   }
   fprintf(cfin,"%s","\n"); 
  }

  fclose(cfin);

  // write the -surf.plt file
  cfin = fopen("cfin-surf.plt","w");

  fprintf(cfin,"%s","TITLE      =  Boundary data\n");
  fprintf(cfin,"%s","VARIABLES  =  \"x0\" \"x1\" ");
  for(unsigned isp=0; isp<(*nsp);isp++) {
   fprintf(cfin,"%s %s",rhostr.at(isp).c_str(), " ");
  }
  fprintf(cfin,"%s","\"u\" \"v\" \"T\" \"Tv0\" \n");

  for(int IBC=0; IBC<maxNCl; IBC++) {

   if(ICLR->at(IBC+1)>0) {

   nbNodestr.str(string());
   Tstr.str(string());
   nbElemstr.str(string());

    ++BND;
    if((IBC+1)==10) {
     if (m_boundary == "single") {
      Tstr << "T= \"" << 10 << ", TR 0\", ";
      nbElemstr << "E= " << ICLR->at(IBC+1) << ", ";
      elemVector.resize(ICLR->at(IBC+1)*2);
      elemVector_noduplicate.resize(ICLR->at(IBC+1)*2);
      elemVector_unitIndex.resize(ICLR->at(IBC+1)*2);
      h=0;
      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         elemVector.at(h) = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1);
         h++;
        }
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

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

      elemVector_sort=elemVector_noduplicate;

      sort(elemVector_sort.begin(),elemVector_sort.end());

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       array_elemVector_sort(0,j)=elemVector_sort.at(j);
       array_elemVector_sort(1,j)=j+1;
      }

      for(unsigned j=0;j<elemVector.size(); j++) {
       for(unsigned k=0;k<array_elemVector_sort.getnCols();k++) {
        if(elemVector.at(j)==array_elemVector_sort(0,k)) {
         elemVector_unitIndex.at(j)=array_elemVector_sort(1,k); }
       }
      }

      nbNodestr << "N= " << elemVector_noduplicate.size() << ", ";

      fprintf(cfin,"%s %s","ZONE",nbNodestr.str().c_str());
      fprintf(cfin,"%s %s",Tstr.str().c_str(),nbElemstr.str().c_str());
      fprintf(cfin,"%s","F=FEPOINT, ET=LINESEG, SOLUTIONTIME=0\n");

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       for(unsigned k=0;k<PhysicsInfo::getnbDim();k++) {
        fprintf(cfin,"%20.16F %s",(*XY)(k,elemVector_sort.at(j)-1)," ");
       }
       for(unsigned k=0;k<(*ndof);k++) {
        fprintf(cfin,"%20.16F %s",(*zroe)(k,elemVector_sort.at(j)-1)," ");
       }
       fprintf(cfin,"%s","\n");
      }
      
      for(unsigned j=0;j<elemVector.size()-1;j=j+2) {
       fprintf(cfin,"%11i",elemVector_unitIndex.at(j));
       fprintf(cfin,"%11i %s",elemVector_unitIndex.at(j+1),"\n");
      }

     } // if m_boundary = single

     if (m_boundary == "splitted") {
      Tstr << "T= InnerSup, TR 0\", ";
      nbElemstr << "E= " << nbSh/2-1 << ", ";
      elemVector.resize((nbSh/2-1)*2);
      elemVector_noduplicate.resize((nbSh/2-1)*2);
      elemVector_unitIndex.resize((nbSh/2-1)*2);

      h=0;
      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         np.at(k) = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1)-1;
        }
        if ((np.at(0) >= minSh) && (np.at(0) <  (minSh+nbSh/2)) &&
            (np.at(1) >= minSh) && (np.at(1) <  (minSh+nbSh/2))) {
         elemVector.at(h)=np.at(0);
         elemVector.at(h+1)=np.at(1);
         h=h+2;
        } // if np conditions
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

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

      elemVector_sort=elemVector_noduplicate;

      sort(elemVector_sort.begin(),elemVector_sort.end());

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       array_elemVector_sort(0,j)=elemVector_sort.at(j);
       array_elemVector_sort(1,j)=j+1;
      }

      for(unsigned j=0;j<elemVector.size(); j++) {
       for(unsigned k=0;k<array_elemVector_sort.getnCols();k++) {
        if(elemVector.at(j)==array_elemVector_sort(0,k)) {
         elemVector_unitIndex.at(j)=array_elemVector_sort(1,k); }
       }
      }

      nbNodestr << "N= " << elemVector_noduplicate.size() << ", ";

      fprintf(cfin,"%s %s","ZONE",nbNodestr.str().c_str());
      fprintf(cfin,"%s %s",Tstr.str().c_str(),nbElemstr.str().c_str());
      fprintf(cfin,"%s","F=FEPOINT, ET=LINESEG, SOLUTIONTIME=0\n");

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       for(unsigned k=0;k<PhysicsInfo::getnbDim();k++) {
        fprintf(cfin,"%20.16F %s",(*XY)(k,elemVector_sort.at(j)-1)," ");
       }
       for(unsigned k=0;k<(*ndof);k++) {
        fprintf(cfin,"%20.16F %s",(*zroe)(k,elemVector_sort.at(j)-1)," ");
       }
       fprintf(cfin,"%s","\n");
      }

      for(unsigned j=0;j<elemVector.size()-1;j=j+2) {
       fprintf(cfin,"%11i",elemVector_unitIndex.at(j));
       fprintf(cfin,"%11i %s",elemVector_unitIndex.at(j+1),"\n");
      }

      nbNodestr.str(string());
      Tstr.str(string());
      nbElemstr.str(string());

      Tstr << "T= InnerSub, TR 0\", ";
      nbElemstr << "E= " << nbSh/2-1 << ", ";
      elemVector.resize((nbSh/2-1)*2);
      elemVector_noduplicate.resize((nbSh/2-1)*2);
      elemVector_unitIndex.resize((nbSh/2-1)*2);

      h=0;
      for(unsigned j=0; j<nbfac->at(1); j++) {
       if((*bndfac)(2,j)==(IBC+1)) {
        int elem = (*bndfac)(0,j);
        int vert = (*bndfac)(1,j);
        for(unsigned k=0; k<2; k++) {
         np.at(k) = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1)-1;
        }
        if ((np.at(0) > (maxSh-nbSh/2)) && (np.at(0) <= maxSh) &&
            (np.at(1) > (maxSh-nbSh/2)) && (np.at(1) <= maxSh)) {
         elemVector.at(h)=np.at(0);
         elemVector.at(h+1)=np.at(1);
         h=h+2;
        } // if np conditions
       } // if (*bndfac)(2,j)==(IBC+1)
      } // for j<nbfac->at(1)

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

      elemVector_sort=elemVector_noduplicate;

      sort(elemVector_sort.begin(),elemVector_sort.end());

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       array_elemVector_sort(0,j)=elemVector_sort.at(j);
       array_elemVector_sort(1,j)=j+1;
      }

      for(unsigned j=0;j<elemVector.size(); j++) {
       for(unsigned k=0;k<array_elemVector_sort.getnCols();k++) {
        if(elemVector.at(j)==array_elemVector_sort(0,k)) {
         elemVector_unitIndex.at(j)=array_elemVector_sort(1,k); }
       }
      }

      nbNodestr << "N= " << elemVector_noduplicate.size() << ", ";

      fprintf(cfin,"%s %s","ZONE",nbNodestr.str().c_str());
      fprintf(cfin,"%s %s",Tstr.str().c_str(),nbElemstr.str().c_str());
      fprintf(cfin,"%s","F=FEPOINT, ET=LINESEG, SOLUTIONTIME=0\n");

      for(unsigned j=0; j<elemVector_sort.size(); j++) {
       for(unsigned k=0;k<PhysicsInfo::getnbDim();k++) {
        fprintf(cfin,"%20.16F %s",(*XY)(k,elemVector_sort.at(j)-1)," ");
       }
       for(unsigned k=0;k<(*ndof);k++) {
        fprintf(cfin,"%20.16F %s",(*zroe)(k,elemVector_sort.at(j)-1)," ");
       }
       fprintf(cfin,"%s","\n");
      }

      for(unsigned j=0;j<elemVector.size()-1;j=j+2) {
       fprintf(cfin,"%11i",elemVector_unitIndex.at(j));
       fprintf(cfin,"%11i %s",elemVector_unitIndex.at(j+1),"\n");
      }

     } // if m_boundary = splitted
    } // if IBC==10

    else {

     Tstr << "T= \"" << BND << ", TR 0\", ";
     nbElemstr << "E= " << ICLR->at(IBC+1) << ", ";
     elemVector.resize(ICLR->at(IBC+1)*2);
     elemVector_noduplicate.resize(ICLR->at(IBC+1)*2);
     elemVector_unitIndex.resize(ICLR->at(IBC+1)*2);
     h=0;
     for(unsigned j=0; j<nbfac->at(1); j++) {
      if((*bndfac)(2,j)==(IBC+1)) {
       int elem = (*bndfac)(0,j);
       int vert = (*bndfac)(1,j);
       for(unsigned k=0; k<2; k++) {
        elemVector.at(h) = (*celnod)(J.callJcycl(vert+k+1)-1,elem-1);
        h++;
       }
      } // if (*bndfac)(2,j)==(IBC+1)
     } // for j<nbfac->at(1)

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

     elemVector_sort=elemVector_noduplicate;

     sort(elemVector_sort.begin(),elemVector_sort.end());

     for(unsigned j=0; j<elemVector_sort.size(); j++) {
      array_elemVector_sort(0,j)=elemVector_sort.at(j);
      array_elemVector_sort(1,j)=j+1;
     }

     for(unsigned j=0;j<elemVector.size(); j++) {
      for(unsigned k=0;k<array_elemVector_sort.getnCols();k++) {
       if(elemVector.at(j)==array_elemVector_sort(0,k)) {
        elemVector_unitIndex.at(j)=array_elemVector_sort(1,k); }
      }
     }

     nbNodestr << "N= " << elemVector_noduplicate.size() << ", ";

     fprintf(cfin,"%s %s","ZONE",nbNodestr.str().c_str());
     fprintf(cfin,"%s %s",Tstr.str().c_str(),nbElemstr.str().c_str());
     fprintf(cfin,"%s","F=FEPOINT, ET=LINESEG, SOLUTIONTIME=0\n");

     for(unsigned j=0; j<elemVector_sort.size(); j++) {
      for(unsigned k=0;k<PhysicsInfo::getnbDim();k++) {
       fprintf(cfin,"%20.16F %s",(*XY)(k,elemVector_sort.at(j)-1)," ");
      }
      for(unsigned k=0;k<(*ndof);k++) {
       fprintf(cfin,"%20.16F %s",(*zroe)(k,elemVector_sort.at(j)-1)," ");
      }
      fprintf(cfin,"%s","\n");
     }

     for(unsigned j=0;j<elemVector.size()-1;j=j+2) {
      fprintf(cfin,"%11i",elemVector_unitIndex.at(j));
      fprintf(cfin,"%11i %s",elemVector_unitIndex.at(j+1),"\n");
     }
    } // else (ICLR(IBC+1)!=10)

   } // if ICLR.at(IBC+1)>0
  } // for IBC<maxNCl

  fclose(cfin);
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::freeArray()
{
  delete XY; delete zroe;
  delete celnod; delete celcel; delete bndfac;
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  celcelVect = MeshData::getInstance().getData <vector<int> >("CELCEL");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  fname = MeshData::getInstance().getData <stringstream>("FNAME");
  ICLR = MeshData::getInstance().getData <vector<int> >("ICLR");
}

//----------------------------------------------------------------------------//

void Triangle2Tecplot::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp  = PhysicsData::getInstance().getData <unsigned> ("NSP");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
