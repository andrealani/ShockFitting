// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/CoPntDispl.hh"
#include "RemeshingSF/CoIntrPnt.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CoPntDispl, Remeshing> coPointdisplProv("CoPntDispl");

//--------------------------------------------------------------------------//

CoPntDispl::CoPntDispl (const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

CoPntDispl::~CoPntDispl()
{
}

//--------------------------------------------------------------------------//

void CoPntDispl::setup()
{
  LogToScreen(VERBOSE, "CoPntDispl::setup() => start\n");

  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "CoPntDispl::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CoPntDispl::unsetup()
{
  logfile.Close();
  LogToScreen(VERBOSE, "CoPntDispl::unsetup()");
}

//--------------------------------------------------------------------------//

void CoPntDispl::remesh()
{
  LogToScreen(INFO, "CoPntDispl::remesh()\n");

  setMeshData();
  setPhysicsData();
  setAddress();

  // initialize nodcodsh values to -99
  setNodCodSh();

  // add second layer of shock nodes
  // place extra shock nodes on the shock normal
  addSecondShockLayer();

  // correct the displacement in the boundary special points
  for(unsigned ISPPNTS=0; ISPPNTS<(*r_nSpecPoints); ISPPNTS++) {

   if (r_typeSpecPoints->at(ISPPNTS) == "IPX" ||
       r_typeSpecPoints->at(ISPPNTS) == "OPX" ||
       r_typeSpecPoints->at(ISPPNTS) == "WPNRX")
    { setCoorForIPXorOPXorWPNRX(ISPPNTS); } 

   else if (r_typeSpecPoints->at(ISPPNTS) == "IPY" ||
            r_typeSpecPoints->at(ISPPNTS) == "OPY" ||
            r_typeSpecPoints->at(ISPPNTS) == "WPNRY")
    { setCoorForIPYorOPYorWPNRY(ISPPNTS); }

   else if (r_typeSpecPoints->at(ISPPNTS) == "TP")
    { setCoorForTP(ISPPNTS); }

   else if (r_typeSpecPoints->at(ISPPNTS) == "RRX")
    { setCoorForRRX(ISPPNTS); }

   else if (r_typeSpecPoints->at(ISPPNTS) == "QP")
    { setCoorForQP(ISPPNTS); }

   else if (r_typeSpecPoints->at(ISPPNTS) == "EP")
    { setCoorForEP(ISPPNTS); }

   else if (r_typeSpecPoints->at(ISPPNTS) == "C")
    { setCoorForC(ISPPNTS); }

   else 
    { logfile ("Condition not implemented\n");
      cout << "Condition not implemented\n";
      exit(1); }
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl::setNodCodSh()
{
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    (*r_NodCodSh)(I,ISH) = -99;
   }
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl::addSecondShockLayer()
{
  for(unsigned ish=0; ish<(*r_nShocks); ish++) {
   for(unsigned i=0; i<r_nShockPoints->at(ish); i++) {
    dx = (*r_vShNor)(0,i,ish);
    dy = (*r_vShNor)(1,i,ish);
    (*r_XYShu)(0,i,ish) = (*r_XYSh)(0,i,ish) + 0.5 * (*eps) * dx;
    (*r_XYShu)(1,i,ish) = (*r_XYSh)(1,i,ish) + 0.5 * (*eps) * dy;
    (*r_XYShd)(0,i,ish) = (*r_XYSh)(0,i,ish) - 0.5 * (*eps) * dx;
    (*r_XYShd)(1,i,ish) = (*r_XYSh)(1,i,ish) - 0.5 * (*eps) * dy;

    (*r_NodCodSh)(i,ish) = 10;
   }
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForIPXorOPXorWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  dx = (*r_vShNor)(0,IP.at(0),ISH.at(0));
  dy = (*r_vShNor)(1,IP.at(0),ISH.at(0));

  (*r_XYShu)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0))
                                      + 0.5 * (*eps)/dx;
  (*r_XYShu)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0));
  (*r_XYShd)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0))
                                      - 0.5 * (*eps)/dx;
  (*r_XYShd)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForIPYorOPYorWPNRY(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  dx = (*r_vShNor)(0,IP.at(0),ISH.at(0));
  dy = (*r_vShNor)(1,IP.at(0),ISH.at(0));

  (*r_XYShu)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0));
  (*r_XYShu)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0)) 
                                     + 0.5 * (*eps)/dy;
  (*r_XYShd)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0));
  (*r_XYShd)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0)) 
                                     - 0.5 * (*eps)/dy;
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForTP(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  // ISH.at(2) mach stem
  // ISH.at(3) contact discontinuity
  setShockIndeces(4,ISPPNTS);

  // move incident shock
  moveDiscontinuity(IP.at(0),ISH.at(0));

  // move reflected shock
  moveDiscontinuity(IP.at(1),ISH.at(1));

  // move mach stem
  moveDiscontinuity(IP.at(2),ISH.at(2));

  // move contact discontinuity
  moveDiscontinuity(IP.at(3),ISH.at(3));

  // superimpose the upstream point of the incident shock to
  // the upstream point of the mach stem
  superimposeDiscPoints("up",IP.at(0),ISH.at(0),
                        "up",IP.at(2),ISH.at(2));

  // superimpose the downstream point of the incident shock to
  // the upstream point of the reflected shock
  superimposeDiscPoints("down",IP.at(0),ISH.at(0),
			"up",  IP.at(1),ISH.at(1));

  // superimpose the downstream point of the mach stem to
  // the downstream point of the contact discontinuity
  // (pay attention to the c.d. direction)
  superimposeDiscPoints("down",IP.at(2),ISH.at(2),
                        "down",IP.at(3),ISH.at(3));

  // superimpose downstream point of the reflected shock to
  // upstream point of the contact discontinuity
  // (pay attention to the c.d. direction)
  superimposeDiscPoints("down",IP.at(1),ISH.at(1),
		        "up",  IP.at(3),ISH.at(3));
}

//--------------------------------------------------------------------------//

void CoPntDispl::superimposeDiscPoints(string discState1, unsigned IP1,
                                       unsigned ISH1,
				       string discState2, unsigned IP2,
                                       unsigned ISH2)
{
  Array3D <double>* XYSh1 = new Array3D <double>(*ndim,*npshmax,*nshmax);
  Array3D <double>* XYSh2 = new Array3D <double>(*ndim,*npshmax,*nshmax);

  if (discState1=="up")   {XYSh1 = r_XYShu;}
  if (discState1=="down") {XYSh1 = r_XYShd;}
  if (discState2=="up")   {XYSh2 = r_XYShu;}
  if (discState2=="down") {XYSh2 = r_XYShd;}

  for(unsigned K=0; K<(*ndim); K++) {
   dum = 0.5 * ((*XYSh1)(K,IP1,ISH1)+(*XYSh2)(K,IP2,ISH2));
   (*XYSh1)(K,IP1,ISH1) = dum;
   (*XYSh2)(K,IP2,ISH2) = dum;
  } 
  delete XYSh1; delete XYSh2;
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForRRX(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  for (unsigned i=0; i<2; i++) {
   (*r_XYShu)(0,IP.at(i),ISH.at(i)) = (*r_XYSh)(0,IP.at(i),ISH.at(i)) 
				      + 0.5 * (*eps)/dx;
   (*r_XYShu)(1,IP.at(i),ISH.at(i)) = (*r_XYSh)(1,IP.at(i),ISH.at(i)); 

   (*r_XYShd)(0,IP.at(i),ISH.at(i)) = (*r_XYSh)(0,IP.at(i),ISH.at(i))
                                      - 0.5 * (*eps)/dx;
   (*r_XYShd)(1,IP.at(i),ISH.at(i)) = (*r_XYSh)(1,IP.at(i),ISH.at(i));
  }

  xc = getDownXVect(IP.at(0),ISH.at(0));
  yc = getDownYVect(IP.at(0),ISH.at(0));

  xs = getUpXVect(IP.at(1),ISH.at(1));
  ys = getUpYVect(IP.at(1),ISH.at(1));

  CoIntrPnt findNewCoor;
  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);
  
  (*r_XYShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  (*r_XYShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY();

  (*r_XYShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  (*r_XYShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY();
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForQP(unsigned ISPPNTS)
{
  // ISH.at(0) first incident shock
  // ISH.at(1) first reflected shock
  // ISH.at(2) second incident shock
  // ISH.at(3) second reflected shock
  // ISH.at(4) contact discontinuity
  setShockIndeces(5,ISPPNTS);

  // set the family which the first incident shock belongs to
  f1 = (*r_vShNor)(0,IP.at(0),ISH.at(0))* (*r_ZRoeShu)(3,IP.at(0),ISH.at(0))
      -(*r_vShNor)(1,IP.at(0),ISH.at(0))* (*r_ZRoeShu)(2,IP.at(0),ISH.at(0)); 
  f1 = copysign(1,f1);

  // set the family which the second incident shock belongs to
  f3 = (*r_vShNor)(0,IP.at(2),ISH.at(2))* (*r_ZRoeShu)(3,IP.at(2),ISH.at(2))
      -(*r_vShNor)(1,IP.at(2),ISH.at(2))* (*r_ZRoeShu)(2,IP.at(2),ISH.at(2));
  f3 = copysign(1,f1);

  // move point belongs to the first incident shock
  moveDiscontinuity(IP.at(0),ISH.at(0));

  // move point belongs to the first incident shock
  moveDiscontinuity(IP.at(1),ISH.at(1));

  // move point belongs to the second incident shock
  moveDiscontinuity(IP.at(2),ISH.at(2));

  // move point belongs to the second reflected shock
  moveDiscontinuity(IP.at(3),ISH.at(3));

  // move point belongs to contact discontinuity
  moveDiscontinuity(IP.at(4),ISH.at(4));

  CoIntrPnt findNewCoor;

  // superimpose point belongs to upstream state of
  // the first incident shock to point belongs to 
  // downstream state of the second incident shock
  // if the incident shocks belong to different families
  if (f1>0) { xc = getUpXVect(IP.at(0),ISH.at(0));
   	      yc = getUpYVect(IP.at(0),ISH.at(0)); }

  else 	    { xc = getDownXVect(IP.at(0),ISH.at(0));
              yc = getDownYVect(IP.at(0),ISH.at(0)); }

  if (f3<0) { xs = getUpXVect(IP.at(2),ISH.at(2));
              ys = getUpYVect(IP.at(2),ISH.at(2)); }

  else      { xs = getDownXVect(IP.at(2),ISH.at(2));
              ys = getDownYVect(IP.at(2),ISH.at(2)); }

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);
  
  if (f1>0) { (*r_XYShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
              (*r_XYShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  else      { (*r_XYShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
              (*r_XYShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  if (f3<0) { (*r_XYShu)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
              (*r_XYShu)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  else      { (*r_XYShd)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
  	      (*r_XYShd)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  // superimpose point belongs to downstream state of first incident
  // shock to point belongs to downstream state of first
  // reflected shock
  if (f1>0) { xc = getDownXVect(IP.at(0),ISH.at(0));
              yc = getDownYVect(IP.at(0),ISH.at(0)); }

  else      { xc = getUpXVect(IP.at(0),ISH.at(0));
   	      yc = getUpYVect(IP.at(0),ISH.at(0)); }
  
  xs = getUpXVect(IP.at(1),ISH.at(1));
  ys = getUpYVect(IP.at(1),ISH.at(1));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  if(f1>0) { (*r_XYShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  	     (*r_XYShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  else 	   { (*r_XYShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
             (*r_XYShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  (*r_XYShu)(0,IP.at(1),ISH.at(1)) = findNewCoor.getX();
  (*r_XYShu)(1,IP.at(1),ISH.at(1)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of second
  // incident shock to point belongs to upstream state of
  // second reflected shock
  if (f3>0) { xc = getUpXVect(IP.at(2),ISH.at(2));
   	      yc = getUpYVect(IP.at(2),ISH.at(2)); }

  else      { xc = getDownXVect(IP.at(2),ISH.at(2));
 	      yc = getDownYVect(IP.at(2),ISH.at(2)); }

  xs = getUpXVect(IP.at(3),ISH.at(3));
  ys = getUpYVect(IP.at(3),ISH.at(3));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  if (f3<0) { (*r_XYShd)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
   	      (*r_XYShd)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  else 	    { (*r_XYShu)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
   	      (*r_XYShu)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  (*r_XYShu)(0,IP.at(3),ISH.at(3)) = findNewCoor.getX();
  (*r_XYShu)(1,IP.at(3),ISH.at(3)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of first
  // reflected shock to point belongs to upstream state of
  // contact discontinuity
  xc = getDownXVect(IP.at(1),ISH.at(1));
  yc = getDownYVect(IP.at(1),ISH.at(1));

  xs = getUpXVect(IP.at(4),ISH.at(4));
  ys = getUpYVect(IP.at(4),ISH.at(4));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  (*r_XYShd)(0,IP.at(1),ISH.at(1)) = findNewCoor.getX();
  (*r_XYShd)(1,IP.at(1),ISH.at(1)) = findNewCoor.getY();

  (*r_XYShu)(0,IP.at(4),ISH.at(4)) = findNewCoor.getX();
  (*r_XYShu)(1,IP.at(4),ISH.at(4)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of 
  // second reflected shock to point belong to downstream
  // state of contact discontinuity
  xc = getDownXVect(IP.at(3),ISH.at(3));
  yc = getDownYVect(IP.at(3),ISH.at(3));

  xs = getDownXVect(IP.at(4),ISH.at(4));
  ys = getDownYVect(IP.at(4),ISH.at(4));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  (*r_XYShd)(0,IP.at(3),ISH.at(3)) = findNewCoor.getX();
  (*r_XYShd)(1,IP.at(4),ISH.at(4)) = findNewCoor.getY();
}

//--------------------------------------------------------------------------/e {

void CoPntDispl::setCoorForEP(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  // superimpose point belongs to downstream state to point
  // belongs to upstream state
  (*r_XYShd)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0));
  (*r_XYShd)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0));
  (*r_XYShu)(0,IP.at(0),ISH.at(0)) = (*r_XYSh)(0,IP.at(0),ISH.at(0));
  (*r_XYShu)(1,IP.at(0),ISH.at(0)) = (*r_XYSh)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

void CoPntDispl::setCoorForC(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  // superimpose point belongs to downstream state to point
  // belongs to upstream state of connection point
  (*r_XYShd)(0,IP.at(1),ISH.at(1)) = (*r_XYShd)(0,IP.at(0),ISH.at(0));
  (*r_XYShd)(1,IP.at(1),ISH.at(1)) = (*r_XYShd)(1,IP.at(0),ISH.at(0));
  (*r_XYShu)(0,IP.at(1),ISH.at(1)) = (*r_XYShu)(0,IP.at(0),ISH.at(0));
  (*r_XYShu)(1,IP.at(1),ISH.at(1)) = (*r_XYShu)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl::getUpXVect(unsigned IP, unsigned ISH)
{
  std::vector <double> x(2);
  x.at(0) = (*r_XYShu)(0,IP,ISH);
  if (IP==0) { x.at(1) = (*r_XYShu)(0,IP+1,ISH); }
  else       { x.at(1) = (*r_XYShu)(0,IP-1,ISH); }
  return x;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl::getUpYVect(unsigned IP, unsigned ISH)
{
  std::vector <double> y(2);
  y.at(0) = (*r_XYShu)(1,IP,ISH);
  if (IP==0) { y.at(1) = (*r_XYShu)(1,IP+1,ISH); }
  else       { y.at(1) = (*r_XYShu)(1,IP-1,ISH); }
  return y;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl::getDownXVect(unsigned IP, unsigned ISH)
{
  std::vector <double> x(2);
  x.at(0) = (*r_XYShd)(0,IP,ISH);
  if (IP==0) { x.at(1) = (*r_XYShd)(0,IP+1,ISH); }
  else       { x.at(1) = (*r_XYShd)(0,IP-1,ISH); }
  return x;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl::getDownYVect(unsigned IP, unsigned ISH)
{ 
  std::vector <double> y(2);
  y.at(0) = (*r_XYShd)(1,IP,ISH);
  if (IP==0) { y.at(1) = (*r_XYShd)(1,IP+1,ISH); }
  else       { y.at(1) = (*r_XYShd)(1,IP-1,ISH); }
  return y;
}

//--------------------------------------------------------------------------//

void CoPntDispl::moveDiscontinuity(unsigned IP, unsigned ISH)
{
  if(IP==1) {
   tx = (*r_XYSh)(0,IP,ISH) - (*r_XYSh)(0,IP-1,ISH);
   ty = (*r_XYSh)(1,IP,ISH) - (*r_XYSh)(1,IP-1,ISH);
  }
  else {
   tx = (*r_XYSh)(0,IP-1,ISH) - (*r_XYSh)(0,IP,ISH);
   ty = (*r_XYSh)(1,IP-1,ISH) - (*r_XYSh)(1,IP,ISH);
  }
  dum = sqrt(pow(tx,2)+pow(ty,2));
  tx = tx/dum;
  ty = ty/dum;

  (*r_XYShu)(0,IP,ISH) = (*r_XYShu)(0,IP,ISH) + (*eps) * tx;
  (*r_XYShu)(1,IP,ISH) = (*r_XYShu)(1,IP,ISH) + (*eps) * ty;

  (*r_XYShd)(0,IP,ISH) = (*r_XYShd)(0,IP,ISH) + (*eps) * tx;
  (*r_XYShd)(1,IP,ISH) = (*r_XYShd)(1,IP,ISH) + (*eps) * ty;
}

//--------------------------------------------------------------------------//

void CoPntDispl::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
{
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*r_SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*r_SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (r_nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl::setAddress()
{
  unsigned start;
  start = (*npoin);
  r_NodCodSh = new Array2D <int> ((*npshmax),(*nshmax),&nodcod->at(start));
  start = (*npoin)*(*ndof);
  r_ZRoeShu =
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
  start = (*npoin) * (*ndim);
  r_XYShu =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
  start = (*npoin) * (*ndim) + (*npshmax) * (*nshmax) * (*ndim);
  r_XYShd =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
}

//--------------------------------------------------------------------------//

void CoPntDispl::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  eps = MeshData::getInstance().getData <double> ("EPS");
  nodcod = MeshData::getInstance().getData <vector<int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coor = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void CoPntDispl::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  r_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  r_nSpecPoints = 
        PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  r_nShockPoints =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  r_nShockEdges =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  r_typeSh =
        PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  r_typeSpecPoints =
        PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  r_XYSh =
        PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
  r_SHinSPPs =
        PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
  r_vShNor =
        PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting


