// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "MeshGeneratorSF/Tricall.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Jcycl.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//
extern void triangulate(char *triswitches, struct triangulateio *in,
			struct triangulateio *out, struct triangulateio *voronoi);

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Tricall, MeshGenerator> triangleCallProv("Tricall");

//--------------------------------------------------------------------------//

Tricall::Tricall(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

Tricall::~Tricall()
{
}

//--------------------------------------------------------------------------//

void Tricall::setup()
{
  LogToScreen(VERBOSE, "Tricall::setup() => start\n");

  LogToScreen(VERBOSE, "Tricall::setup() => end\n");
}

//--------------------------------------------------------------------------//

void Tricall::unsetup()
{
  LogToScreen(VERBOSE, "Tricall::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Tricall::generate()
{
  LogToScreen(INFO, "Tricall::generate()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  int in1, in2;

  // set ICLR vector to 0 value
  ICLR->assign(20,0);

  triangleReport = fopen("./log/TriangleReport.log","w");

  ilist = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                             PhysicsInfo::getnbShPointsMax();

  // M02M1 and M12M0 are filled with indeces that start 
  // from 1 to NPOIN+2*NSHMAX*NPSHMAX+1
  M02M1->resize(ilist+1); // c++ indeces start from 0
  M12M0->resize(ilist+1); // c++ indeces start from 0

  TNPOIN = 0;

  // set map vector for nodcod
  setMapVectorForNodcod();

  icount = npoin->at(0);

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  ilist = TNPOIN;

  in.numberofpoints = ilist;
  in.numberofpointattributes = (*ndof);

  // assign the size to the triangulation vectors
  setTriSize();

  // store mesh points coordinates and state for the triangulation
  storeMeshVariables();

  icount = npoin->at(0);

  // store upstream shock points coordinates and states for the triangulation
  storeUpstreamStatus();

  // store downstream shock points coordinates and states for the triangulation
  storeDownstreamStatus();

  // set map vector for bndfac and store bndfac for the triangulation
  storeBndfac();

  // compute number of holes
  computenbHoles();

  fprintf(triangleReport,"%s","Input point set:\n\n");

  report(&in, 1, 0, 0, 1, 0, 0);

  // Make necessary initializations so that Triangle can return
  // a triangulation in `mid'
  midInitialization();

  // Triangulate the points.  Switches are chosen to read and write a 
  // PSLG (p), produce an edge list (e) and a triangle neighbor list (n).
  // Switch Q is chosen to suppress all explanation of what Triangle is doing,
  // unless an error occurs
  triangulate((char *)("nepQ"),&in,&mid,(struct triangulateio *) NULL);

  fprintf(triangleReport,"Initial triangulation:\n\n");

  report(&mid, 1, 1, 1, 1, 1, 0);

  // Make necessary initializations so that Triangle can return
  // a triangulation in `out'.
  outInitialization();                                              

  // Refine the triangulation
  // Switch (r) refines a previously generated mesh
  triangulate((char*)("neprQ"), &mid, &out, (struct triangulateio *) NULL); 

  fprintf(triangleReport,"Refined triangulation:\n\n");

  report(&out, 1, 1, 1, 1, 1, 0);

  // Assign the new values computed by triangle

  npoin->at(1) = out.numberofpoints;
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
  // assign mesh points coordinates and state
  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   unsigned i = IPOIN;
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
    (*XY)(IA,IPOIN) = out.pointlist[i*PhysicsInfo::getnbDim()+IA]; }
   unsigned j = IPOIN;
   for(unsigned IA=0; IA<(*ndof); IA++) {
    (*zroe)(IA,IPOIN) = out.pointattributelist[j*(*ndof)+IA]; }
  }

  nelem->at(1) = out.numberoftriangles;
  // resize arrays whit the new nelem value
  totsize = nelem->at(0) + nelem->at(1);
  celcelVect->resize((*nvt) * totsize);
  celnodVect->resize((*nvt) * totsize);
  start = 3 * nelem->at(0);
  celnod = new Array2D<int> (3, nelem->at(1), &celnodVect->at(start));
  celcel = new Array2D<int> (3, nelem->at(1), &celcelVect->at(start));
  // assign celnod array values
  for(unsigned IELEM=0; IELEM<nelem->at(1); IELEM++) {
   unsigned i = IELEM;
   for(unsigned J=0; J<(*nvt); J++) {
    (*celnod)(J,IELEM) = out.trianglelist[i*(*nvt)+J]; 
    (*celcel)(J,IELEM) = out.neighborlist[i*(*nvt)+J]; 
   }
  }

  nbfac->at(1) = out.numberofsegments;
  // resize array with the new nbfac value
  totsize = nbfac->at(0) + nbfac->at(1) +
            4 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax();
  bndfacVect->resize(3 * totsize);
  start = 3 * (nbfac->at(0) +
          2 * PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax());
  bndfac = new Array2D<int> (3,(nbfac->at(1) +
                             2 * PhysicsInfo::getnbShMax() *
                             PhysicsInfo::getnbShEdgesMax()),
                             &bndfacVect->at(start));

  // fill bndfac values
  // only for element bndfac(0,*) and bndfac(1,*)
  // bndfac(2,*) will be filled after

  // create Jcycl object
  Jcycl J;
  unsigned ibfac=0;
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

  unsigned nbBoundaryfaces = ibfac;

  unsigned nedge = out.numberofedges;

  // set ICLR vector to 0 value
  ICLR->assign(20,0);

  // dummy variables
  int n1, n2, ielem, IBC, ivert;
  unsigned IFACE=0; unsigned I=0; 
  ibfac = 0;

  while(IFACE<nedge) {
   n1 = out.edgelist[IFACE*((*nvt)-1)];
   n2 = out.edgelist[IFACE*((*nvt)-1)+1];
   IBC = out.edgemarkerlist[IFACE];
   // if a boundary faces, look for the parent element in bndfac
   if(IBC==0) { ICLR->at(IBC) = ICLR->at(IBC)+1; }
   else {
    ++ibfac;
    I=0;
    while(I<(nbBoundaryfaces)) {

     ielem = (*bndfac)(0,I);
     ivert = (*bndfac)(1,I);

     in1 = (*celnod)(J.callJcycl(ivert+1)-1,ielem-1);//c++ indeces start from 0
     in2 = (*celnod)(J.callJcycl(ivert+2)-1,ielem-1);//c++ indeces start from 0

     if ((in1==n1 && in2==n2) || (in1==n2 && in2==n1)) {
     if((*bndfac)(2,I)!=-1) { cout << "Tricall::error => " << I;
                              cout << " " << (*bndfac)(2,I) << "\n";
                              exit(1);                                  }
      (*bndfac)(2,I) = IBC;
      ICLR->at(IBC) = ICLR->at(IBC) + 1;
      goto nine;
     } // if (in1==n1 && in2==n2) || (in1==n2 && in2==n1)
     I++;
    } // while I<(nbBoundaryfaces)
    cout << "Tricall::error => Cannot match boundary face ";
    cout << IFACE << " " << n1 << " " << n2 << " " << IBC << "\n";
    exit(IFACE);
    }
nine:
   IFACE++;
   }

  if (ibfac!=nbBoundaryfaces) {
   cout << "Tricall::error => No matching boundary faces\n";
   exit(1);
  }

  // Free all allocated arrays, including those allocated by Triangle.
  freeTri();

  fclose(triangleReport);
}

//--------------------------------------------------------------------------//

void Tricall::generate(string dummystr)
{
  LogToScreen(INFO, "Tricall::generate()\n");




}

//--------------------------------------------------------------------------//

void Tricall::setMapVectorForNodcod()
{
  for (unsigned IPOIN=0; IPOIN< npoin->at(0); IPOIN++) {
   if (nodcod->at(IPOIN)>=0) { ++TNPOIN;
                               M02M1->at(IPOIN+1) = TNPOIN;
                               M12M0->at(TNPOIN) = IPOIN+1; } 
  }
}

//--------------------------------------------------------------------------//

void Tricall::setMapVectorForNodcodSh()
{
  for (unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
    ++icount;

    if ((*NodCodSh)(I,ISH)==10) { ++TNPOIN;
                                    M02M1->at(icount) = TNPOIN;
                                    M12M0->at(TNPOIN) = icount;}
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::storeMeshVariables()
{
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if(nodcod->at(IPOIN)>=0) {
    unsigned i=M02M1->at(IPOIN+1)-1; // triangle vectors are stored from index 0
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     in.pointlist[i*PhysicsInfo::getnbDim()+IA] = (*XY)(IA,IPOIN);}
    unsigned j=M02M1->at(IPOIN+1)-1; // triangle vectors are stored from index 0
    for(unsigned IA=0; IA<(*ndof); IA++) {
     in.pointattributelist[j*(*ndof)+IA] = (*zroe)(IA,IPOIN);}
    in.pointmarkerlist[i] = nodcod->at(IPOIN);
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::storeUpstreamStatus()
{
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
    ++icount;
    if((*NodCodSh)(I,ISH)==10) {
     unsigned i=M02M1->at(icount)-1; // triangle vectors are stored from index 0
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
      in.pointlist[i*PhysicsInfo::getnbDim()+IA] = (*XYShu)(IA,I,ISH); }
     unsigned j=M02M1->at(icount)-1; // triangle vectors are stored from index 0
     for(unsigned IA=0; IA<(*ndof); IA++) {
      in.pointattributelist[j*(*ndof)+IA] = (*ZRoeShu)(IA,I,ISH); }
     // triangle vectors are stored from index 0
     in.pointmarkerlist[i] = (*NodCodSh)(I,ISH);
    }
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::storeDownstreamStatus()
{
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<PhysicsInfo::getnbShPointsMax(); I++) {
     icount++;
    if((*NodCodSh)(I,ISH)==10) {
     unsigned i=M02M1->at(icount)-1; // triangle vectors are stored from index 0
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
      in.pointlist[i*PhysicsInfo::getnbDim()+IA] = (*XYShd)(IA,I,ISH); }
     unsigned j=M02M1->at(icount)-1; // triangle vectors are stored from index 0
     for(unsigned IA=0; IA<(*ndof); IA++) {
      in.pointattributelist[j*(*ndof)+IA] = (*ZRoeShd)(IA,I,ISH); }
     // triangle vectors are stored from index 0
     in.pointmarkerlist[i] = (*NodCodSh)(I,ISH);
    }
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::storeBndfac()
{
  int IBC;

  ICHECK = 0;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(2,IFACE);
   if(IBC>0) { ++ICHECK; }
  }
  in.numberofsegments = ICHECK;
  in.numberofregions = 0;
  ICHECK=0;

  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) { 
   IBC = (*bndfac)(2,IFACE);

   if(IBC>0) {
    ++ICHECK; 
    unsigned i=ICHECK-1; // triangle vectors are stored from index 0
    for(unsigned IA=0; IA<2; IA++) {
     in.segmentlist[i*2+IA] = M02M1->at((*bndfac)(IA,IFACE)); }
    // triangle vectors are stored from index 0
    in.segmentmarkerlist[i] = IBC;
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::computenbHoles()
{
  nHoles = 0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   nHoles = nHoles + nShockPoints->at(ISH)-2;
  }

  in.numberofholes = nHoles + MeshData::getInstance().getnbAddHoles();

  unsigned iHole=0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned I=1; I<nShockPoints->at(ISH)-1; I++) {
    ++iHole;
    unsigned i=iHole-1; // triangle vectors are stored from index 0
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
    in.holelist[i*PhysicsInfo::getnbDim()+IA] = (*XYSh)(IA,I,ISH);
    }
   }
  }

  for (unsigned I=(2*iHole)-1; I<MeshData::getInstance().getnbAddHoles(); I++) {
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
   in.holelist[I*PhysicsInfo::getnbDim()+IA] = caddholes->at(IA);
   }
  }
}

//--------------------------------------------------------------------------//

void Tricall::midInitialization()
{
  mid.pointlist = (double *) NULL; 
  mid.pointattributelist = (double *) NULL;
  mid.pointmarkerlist = (int *) NULL;
  mid.trianglelist = (int *) NULL;  
  mid.triangleattributelist = (double *) NULL;
  mid.neighborlist = (int *) NULL; 
  mid.segmentlist = (int *) NULL;
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;  
  mid.edgemarkerlist = (int *) NULL; 
}

//--------------------------------------------------------------------------//

void Tricall::outInitialization()
{
  out.pointlist = (double *) NULL;
  out.pointattributelist = (double *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL;
  out.neighborlist = (int *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int*) NULL;
  out.edgelist = (int *) NULL;
  out.edgemarkerlist = (int *) NULL;
}

//--------------------------------------------------------------------------//

void Tricall::freeTri()
{
  free(in.pointlist);
  free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(in.holelist);
  free(mid.pointlist);
  free(mid.pointattributelist);
  free(mid.pointmarkerlist);
  free(mid.trianglelist);
  free(mid.triangleattributelist);
  free(mid.neighborlist);
  free(mid.segmentlist);
  free(mid.segmentmarkerlist);
  free(mid.edgelist);
  free(mid.edgemarkerlist);
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.pointmarkerlist);
  free(out.trianglelist);
  free(out.neighborlist);
  free(out.segmentlist);
  free(out.segmentmarkerlist);
  free(out.edgelist);
  free(out.edgemarkerlist);
}

//--------------------------------------------------------------------------//

void Tricall::report(struct triangulateio* io, int markers, 
                     int reporttriangles, int reportneighbors,
                     int reportsegments, int reportedges, int reportnorms)
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    fprintf(triangleReport,"%s %i","Points ", i);
    for (j = 0; j < 2; j++) {
      fprintf(triangleReport,"%s %.16f","  ", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      fprintf(triangleReport,"%s","  attributes ");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      fprintf(triangleReport,"%s %.16f","  ",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      fprintf(triangleReport,"%s %i %s","   marker ", io->pointmarkerlist[i],"\n");
    } else {
      fprintf(triangleReport,"%s","\n");
    }
  }
  fprintf(triangleReport,"%s","\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        fprintf(triangleReport,"%s %i","Triangle points: ", i);
        for (j = 0; j < io->numberofcorners; j++) {
          fprintf(triangleReport,"%s %i","  ", 
                  io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          fprintf(triangleReport,"%s","   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          fprintf(triangleReport,"%s %.16f","  ", io->triangleattributelist[i *
                                          io->numberoftriangleattributes + j]);
        }
        fprintf(triangleReport,"%s","\n");
      }
      if (reportneighbors) {
        fprintf(triangleReport,"%s %i","Triangle neighbors: ", i);
        for (j = 0; j < 3; j++) {
          fprintf(triangleReport,"%s %i","  ", io->neighborlist[i * 3 + j]);
        }
        fprintf(triangleReport,"%s","\n");
      }
    }
    fprintf(triangleReport,"%s","\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      fprintf(triangleReport,"%s %i","Segment points:", i);
      for (j = 0; j < 2; j++) {
        fprintf(triangleReport,"%s %i","  ", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        fprintf(triangleReport,"%s %i %s","   marker ", 
                               io->segmentmarkerlist[i],"\n");
      } else {
        fprintf(triangleReport,"%s","\n");
      }
    }
    fprintf(triangleReport,"%s","\n");

    for (i = 0; i < io->numberofholes; i++) {
      fprintf(triangleReport,"%s %i","Hole :", i);
      for (j = 0; j < 2; j++) {
        fprintf(triangleReport,"%s %0.16f","  ", io->holelist[i * 2 + j]);
      }
      fprintf(triangleReport,"%s","\n");
    }
    fprintf(triangleReport,"%s","\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      fprintf(triangleReport,"%s %i","Edge points:", i);
      for (j = 0; j < 2; j++) {
        fprintf(triangleReport,"%s %i","  ", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          fprintf(triangleReport,"%s %.16f","  ", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
       fprintf(triangleReport,"%s %i %s","   marker ", io->edgemarkerlist[i],"\n");
      } else {
        fprintf(triangleReport,"%s","\n");
      }
    }
    fprintf(triangleReport,"%s","\n");
  }
}

//--------------------------------------------------------------------------//

void Tricall::setTriSize()
{
  // it is the vector storing the mesh points coordinates
  in.pointlist = (double *) malloc((PhysicsInfo::getnbDim() *
                                    (npoin->at(0) + 2 *
                                     PhysicsInfo::getnbShMax() *
                                     PhysicsInfo::getnbShPointsMax())) *
                                    sizeof(double));
  // it is the vector storing the mesh points states
  in.pointattributelist = (double *) malloc(((*ndof) *
                                             (npoin->at(0) + 2 *
                                              PhysicsInfo::getnbShMax() *
                                              PhysicsInfo::getnbShPointsMax())) *
                                            sizeof(double));
  // it is the vector storing the nodcod
  in.pointmarkerlist = (int *) malloc((npoin->at(0) +
                                       PhysicsInfo::getnbShMax() *
                                       PhysicsInfo::getnbShPointsMax()) *
                                      sizeof(int));
  // it is the vector storing the following elements of bndfac
  // bndfac(0,0:nbfac) and bndfac(1,0:nbfac)
  in.segmentlist = (int *) malloc((2 * (nbfac->at(0) + 2 *
                                       PhysicsInfo::getnbShMax() *
                                       PhysicsInfo::getnbShPointsMax())) *
                                  sizeof(int));
  // it is the vector storing the following elements of bndfac
  // bndfac(2,0:nbfac)
  in.segmentmarkerlist = (int *) malloc((nbfac->at(0) + 2 *
                                         PhysicsInfo::getnbShMax() *
                                         PhysicsInfo::getnbShPointsMax()) *
                                        sizeof(int));
  // it is the vector storing the hole points coordinates
  in.holelist = (double *) malloc((PhysicsInfo::getnbDim() *
                                   (MeshData::getInstance().getnbAddHoles() +
                                    PhysicsInfo::getnbShMax() *
                                    PhysicsInfo::getnbShPointsMax())) *
                                  sizeof(int));
}

//--------------------------------------------------------------------------//

void Tricall::setAddress()
{
  start = 0;
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),npoin->at(0),
                             &coorVect->at(start));
  zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),npoin->at(0),
                               &zroeVect->at(start));
  totsize = nbfac->at(0) + 2 *
                           PhysicsInfo::getnbShMax() *
                           PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  celnod = new Array2D<int> ((*nvt), nelem->at(0), &celnodVect->at(start));
  celcel = new Array2D<int> ((*nvt), nelem->at(0), &celcelVect->at(start));
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZRoeShu = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() + 
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() ;
  XYShu = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim()  +
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim() ;
  XYShd = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
}

//--------------------------------------------------------------------------//

void Tricall::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem =MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbfacSh = MeshData::getInstance().getData<unsigned>("NBFACSH");
  caddholes = 
    MeshData::getInstance().getData <vector<double> > ("CADDholes");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  celcelVect = MeshData::getInstance().getData <vector<int> >("CELCEL");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
  ICLR = MeshData::getInstance().getData <vector <int> > ("ICLR");
}

//--------------------------------------------------------------------------//

void Tricall::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
