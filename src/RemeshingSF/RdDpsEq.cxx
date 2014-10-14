// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/RdDpsEq.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Remeshing.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<RdDpsEq, Remeshing> rdDpsEqProv("RdDpsEq");

//--------------------------------------------------------------------------//

RdDpsEq::RdDpsEq(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

RdDpsEq::~RdDpsEq()
{
}

//--------------------------------------------------------------------------//

void RdDpsEq::setup()
{
  LogToScreen(VERBOSE, "RdDpsEq::setup() => start\n");

  setPhysicsData();
  setMeshData();
  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "RdDpsEq::setup() => end\n");
}

//--------------------------------------------------------------------------//

void RdDpsEq::unsetup()
{
  LogToScreen(VERBOSE, "RdDpsEq::unsetup()\n");

  logfile.Close();
}

//--------------------------------------------------------------------------//

void RdDpsEq::remesh()
{
  LogToScreen(INFO, "RdDpsEq::remesh()\n");

  // assign start pointers of Array2D and 3D
  setAddress();

  // resize vectors and arrays
  setSize();

  for (unsigned I=0; I<(*r_nShocks); I++){

   // compute length of shock edge
   computeShEdgeLength(I);

   // compute number of shock points of redistribution
   // and check the new number of shock points
   computeNbDistrShPoints(I);

   // compute distribution step
   computeDistrStep(I);

   // interpolation of new shock points
   // set the first and the last point
   newPointsInterp(I);

   // computation of internal points
   computeInterPoints(I);

   // rewrite the state arrays and indices
   rewriteValues(I);
  }
}

//--------------------------------------------------------------------------//

void RdDpsEq::computeShEdgeLength(unsigned index)
{
  ISH=index;
  unsigned ish = ISH+1;
  logfile("Shock n. :", ish,"\n");
  Sh_ABSC.at(0) = 0;
  for (unsigned IV=0; IV<r_nShockEdges->at(ISH); IV++) {
   unsigned I = IV+1;
   Sh_Edge_length = 0;
   for (unsigned K=0; K<(*ndim); K++) {
    Sh_Edge_length = Sh_Edge_length +
                     pow(((*r_XYSh)(K,IV,ISH)-(*r_XYSh)(K,I,ISH)),2);
   }
   Sh_Edge_length = sqrt(Sh_Edge_length);
   Sh_ABSC.at(I) = Sh_ABSC.at(I-1) + Sh_Edge_length;
  }

  logfile("Number of points old distribution:", r_nShockPoints->at(ISH),"\n");
}

//--------------------------------------------------------------------------//

void RdDpsEq::computeNbDistrShPoints(unsigned index)
{
  ISH=index;
  nShockPoints_new = Sh_ABSC.at(r_nShockPoints->at(ISH)-1)/(*r_dxcell)+1;
  if (nShockPoints_new > (*npshmax)) {
   cout << "Too many shock points! Increase NPSHMAX in input.case\n";
   exit(1);}
}

//--------------------------------------------------------------------------//

void RdDpsEq::computeDistrStep(unsigned index)
{
  ISH=index;
  double dx = Sh_ABSC.at(r_nShockPoints->at(ISH)-1)/(nShockPoints_new-1);
  Sh_ABSC_New.at(0) = 0;
  for (unsigned IV=0; IV<nShockPoints_new; IV++) {
   unsigned I = IV+1;
   Sh_ABSC_New.at(I) = Sh_ABSC_New.at(I-1) + dx;

  }
  logfile ("Number of points new distribution:",nShockPoints_new, "\n");
}

//--------------------------------------------------------------------------//

void RdDpsEq::newPointsInterp(unsigned index)
{
  ISH = index;
   for (unsigned K=0; K<(*ndim); K++) {
    XYSh_New(K,0) = (*r_XYSh)(K,0,ISH);
   }
   for (unsigned K=0; K<(*ndof); K++) {
    ZRoeShu_New(K,0) = (*r_ZRoeShu)(K,0,ISH);
    ZRoeShd_New(K,0) = (*r_ZRoeShd)(K,0,ISH);
   }
   for (unsigned K=0; K<(*ndim); K++) {
    XYSh_New(K,nShockPoints_new-1) = (*r_XYSh)(K,r_nShockPoints->at(ISH)-1,ISH);
   }
   for (unsigned K=0; K<(*ndof); K++) {
    ZRoeShu_New(K,nShockPoints_new-1) =
                              (*r_ZRoeShu)(K,r_nShockPoints->at(ISH)-1,ISH);
    ZRoeShd_New(K,nShockPoints_new-1) =
                              (*r_ZRoeShd)(K,r_nShockPoints->at(ISH)-1,ISH);
   }
}

//--------------------------------------------------------------------------//

void RdDpsEq::computeInterPoints(unsigned index)
{
  double ds, dsj, alpha, beta;
  unsigned j1, i1;
  ISH=index;
  for(unsigned i=1; i<nShockPoints_new-1; i++) {
   for (unsigned j=1; j<r_nShockPoints->at(ISH); j++) {
    if ( Sh_ABSC_New.at(i) >= Sh_ABSC.at(j-1) && 
         Sh_ABSC_New.at(i) <  Sh_ABSC.at(j) ) {
     ds = Sh_ABSC.at(j)-Sh_ABSC.at(j-1);
     dsj = Sh_ABSC_New.at(i)-Sh_ABSC.at(j-1);
     alpha = dsj/ds;
     beta = 1-alpha;
     j1 = j-1; i1 = i-1;
     logfile("node=",j1,"\n");
     logfile("node_new=",i1,"\n");
     for (unsigned K=0; K<(*ndim); K++) {
      XYSh_New(K,i) = beta * (*r_XYSh)(K,j1,ISH) + alpha * (*r_XYSh)(K,j,ISH);
     }
     for (unsigned K=0; K<(*ndof); K++) {
      ZRoeShu_New(K,i) = beta * (*r_ZRoeShu)(K,j1,ISH) + 
                         alpha * (*r_ZRoeShu)(K,j,ISH);
      ZRoeShd_New(K,i) = beta * (*r_ZRoeShd)(K,j1,ISH) + 
                         alpha * (*r_ZRoeShd)(K,j,ISH); 
      logfile("ZROESHd_j-1(",K,"=",(*r_ZRoeShd)(K,j1,ISH), "\n");
     } // K
    } // if
   } // j
  } // i
}

//--------------------------------------------------------------------------//

void RdDpsEq::rewriteValues(unsigned index)
{
  ISH=index;
  logfile("new shock point coordinates");
  for (unsigned I=0; I<nShockPoints_new; I++) {
   for(unsigned K=0; K<(*ndim); K++) {
    (*r_XYSh)(K,I,ISH)= XYSh_New(K,I);
    //logfile(K)
    logfile((*r_XYSh)(K,I,ISH), " ");
   }
   logfile("\n");
   for(unsigned IV=0; IV<(*ndof); IV++) {
    logfile(ZRoeShd_New(IV,I), " ");
    (*r_ZRoeShu)(IV,I,ISH) = ZRoeShu_New(IV,I);
    (*r_ZRoeShd)(IV,I,ISH) = ZRoeShd_New(IV,I);
   }
   logfile("\n");
  }
  r_nShockPoints->at(ISH)=nShockPoints_new;
  r_nShockEdges->at(ISH) = nShockPoints_new-1;  
}

//--------------------------------------------------------------------------//

void RdDpsEq::setAddress()
{
  unsigned start;
  start = (*npoin)*(*ndof);
  r_ZRoeShu = new Array3D <double> 
              ((*ndof),(*npshmax),(*nshmax),&r_zroe->at(start));
  start = (*npoin) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  r_ZRoeShd = new Array3D <double> 
              ((*ndofmax),(*npshmax),(*nshmax),&r_zroe->at(start));
}

//--------------------------------------------------------------------------//

void RdDpsEq::setSize()
{
  XYSh_New.resize( (*ndim) , (*npshmax) );
  ZRoeShu_New.resize( (*ndof) , (*npshmax) );
  ZRoeShd_New.resize( (*ndof) , (*npshmax) );
  Sh_ABSC.resize( (*npshmax) );
  Sh_ABSC_New.resize( (*npshmax) );
}

//--------------------------------------------------------------------------//

void RdDpsEq::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  r_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks"); 
  r_nShockPoints = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  r_nShockEdges = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  r_XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

void RdDpsEq::setMeshData()
{
  r_zroe = MeshData::getInstance().getData <vector <double> >("ZROE");
  r_dxcell = MeshData::getInstance().getData <double> ("DXCELL");
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
