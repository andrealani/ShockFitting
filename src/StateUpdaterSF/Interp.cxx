// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/Interp.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Area.hh"
#include "MathTools/Jcycl.hh"
#include "MathTools/MinMax.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Interp, StateUpdater> interpProv("Interp");

//--------------------------------------------------------------------------//

Interp::Interp(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

Interp::~Interp()
{
}

//--------------------------------------------------------------------------//

void Interp::setup()
{
  LogToScreen(VERBOSE,"Interp::setup() => start\n");

  LogToScreen(VERBOSE,"Interp::setup() => end\n");
}

//--------------------------------------------------------------------------//

void Interp::unsetup()
{
  LogToScreen(VERBOSE,"Interp::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Interp::update()
{
  LogToScreen(INFO,"Interp::update()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  unsigned I1, I2, JPOIN;

  logfile.Open(getClassName().c_str());

  logfile("Entering in Interp\n");

  // match upstream nodes with downstream ones
  // these are the coordinates of the shocked mesh (1)
  // the nof of shock points belongs the shocked mesh
  // NOT the background mesh
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned IPOIN=0; IPOIN<nShockPoints->at(ISH); IPOIN++)  {
    I1 = ISH * PhysicsInfo::getnbShMax() + IPOIN; //c++ indeces start from 0
    JPOIN = M02M12->at(I1);
    (*XY)(0,JPOIN-1) = (*XYSh)(0,IPOIN,ISH); // c++ indeces start from 0
    (*XY)(1,JPOIN-1) = (*XYSh)(1,IPOIN,ISH); // c++ indeces start from 0
    I2 = I1 + PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax();
    JPOIN = M02M12->at(I2);
    (*XY)(0,JPOIN-1) = (*XYSh)(0,IPOIN,ISH); // c++ indeces start from 0
    (*XY)(1,JPOIN-1) = (*XYSh)(1,IPOIN,ISH); // c++ indeces start from 0
   }
  }

  // interpolate background grid nodes using shocked grid connectivity
  // the interpolation is necessary only for the phantom nodes
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   unsigned m_IPOIN = IPOIN+1;

   // nodcod = -1 internal phantom nodes
   // nodcod = -2 boundary phantom nodes
   if (   (nodcod->at(IPOIN)==-1) || (nodcod->at(IPOIN)==-2)  ) {  
    logfile("Trying to locate ", m_IPOIN);
    logfile("(",(*XYBkg)(0,IPOIN),",",(*XYBkg)(1,IPOIN),")\n");

    finder(IPOIN);

    logfile("Found in cell ",getCell()," ",getIfound(), "\n");
    if(getIfound()!=0) {
     logfile("Search failed for vertex ", m_IPOIN, "\n");
     for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++)
      { logfile((*XYBkg)(I,IPOIN), " "); }
     logfile ("\ncell no. is ",getCell(),"\n");
     exit(1);
    }
   }
   else {
    // M02M1 is filled with indeces that start from 1
    // M02M1(1:NPOIN+2*NSHMAX*NPSHMAX)
    JPOIN = M02M1->at(m_IPOIN);
    if (JPOIN==0) { 
     cout << "Interp::error => something wrong for " <<  JPOIN <<"\n";
    }

    for(unsigned I=0; I<(*ndof); I++) {
     (*zBkg)(I,IPOIN) = (*zroe)(I,JPOIN-1); // c++ indeces start form 0
    }
   }
  }

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void Interp::finder(unsigned IPOIN)
{
  // x0,y0 are the coordinates of the point to be located
  double x0, y0;
  double s, t, help;
  unsigned I;
  vector<double> xp(4);
  vector<double> yp(4);
  vector<unsigned> idxs(3);
  vector<double> a(3);

  // create object computing area
  Area boundedArea; 

  // create object Jcycl
  Jcycl J;

  // create object finding minimum and maximum values
  MinMax <double> m;

  ifound = 0;

  x0 = (*XYBkg)(0,IPOIN); 
  y0 = (*XYBkg)(1,IPOIN);

  unsigned IELEM;

  for(IELEM=0; IELEM<nelem->at(1); IELEM++) {
   for(unsigned IV=0; IV<3; IV++) {
    I = (*celnod)(IV,IELEM); // global code number
    xp.at(IV) = (*XY)(0,I-1);
    yp.at(IV) = (*XY)(1,I-1);
   }
   xp.at(3) = x0;  // node to be located
   yp.at(3) = y0;  // node to be located

   idxs.at(0) = 0; idxs.at(1) = 1; idxs.at(2) = 2; // c++ indeces start from 0
   boundedArea.computeArea(xp,yp,idxs);
   help = 1/boundedArea.getArea();
   idxs.at(2) = 3; // node to be located

   // compute area coordinates (in the x-y plane)
   for(unsigned IV=0; IV<3; IV++) {
    idxs.at(0) = J.callJcycl(1+IV+1)-1;
    idxs.at(1) = J.callJcycl(2+IV+1)-1;
    boundedArea.computeArea(xp,yp,idxs);
    a.at(IV) = boundedArea.getArea() * help;
   }

   s = m.min(a.at(0), a.at(1), a.at(2));
   t = m.max(a.at(0), a.at(1), a.at(2));

   if (( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)) {

    for(unsigned ivar=0; ivar<(*ndof); ivar++) { (*zBkg)(ivar,IPOIN) = 0; }
    for(unsigned IV=0; IV<3; IV++) {
     I = (*celnod)(IV,IELEM);
     help = a.at(IV);
     for(unsigned ivar=0; ivar<(*ndof); ivar++) {
      (*zBkg)(ivar,IPOIN) = (*zBkg)(ivar,IPOIN) + help * (*zroe)(ivar,I-1);
     }
    }
    ielem = IELEM+1; // c++ indeces start from 0
    ifound = 0;
    return;
   } // if ( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)
  }  // for IELEM<nelem->at(1)

  ifound = 1;
  cout << "Search failed for Vertex coords " << x0 << ", " << y0 << "\n";
  for(unsigned IV=0; IV<3; IV++) { cout << a.at(IV) << " "; }
  cout << a.at(0)+a.at(1)+a.at(2) << " " << s << " " << t << "\n";
  return;
}

//--------------------------------------------------------------------------//

void Interp::setAddress()
{
  unsigned start;
  start = 0;
  XYBkg = new Array2D<double>(PhysicsInfo::getnbDim() , 
                              npoin->at(0) + 2 *
                              PhysicsInfo::getnbShPointsMax() *
                              PhysicsInfo::getnbShMax(), 
                              &coorVect->at(start));
  zBkg = new Array2D<double>(PhysicsInfo::getnbDofMax() ,
                              npoin->at(0) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax(),
                              &zroeVect->at(start));
  start = PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                             &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double>(PhysicsInfo::getnbDofMax(),
                              (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                              &zroeVect->at(start));
  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  // XYShu and XYShd have the starting pointers referred to shocked mesh
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  XYShu = new Array3D <double> (PhysicsInfo::getnbDim(),
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          PhysicsInfo::getnbDim() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  XYShd = new Array3D <double> (PhysicsInfo::getnbDim(),
                                PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  // M02M1 is filled with the indeces that start from 1
  // M02M1(1:NPOIN(0)+2*NSHMAX*NPSHMAX)
  // M02M1 is filled wth the indeces that start from 0
  // M02M12(0:2*NSHMAX*NPSHMAX-1)
  M02M12 = new vector<int>(2 * PhysicsInfo::getnbShMax() *
                               PhysicsInfo::getnbShPointsMax());

  for(unsigned i=0; i<M02M12->size(); i++) {
   M02M12->at(i) = M02M1->at(i+npoin->at(0)+1);
  }
}

//--------------------------------------------------------------------------//

void Interp::freeArray()
{
  delete zroe; delete XY;
  delete XYShu; delete XYShd;
  delete celnod;
  delete XYBkg; delete zBkg;
  delete M02M12;
}

//--------------------------------------------------------------------------//

void Interp::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  zroeVect = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
  celnodVect = 
    MeshData::getInstance().getData <vector<int> > ("CELNOD");
  nodcod = 
    MeshData::getInstance().getData <vector<int> > ("NODCOD");
  M02M1 = MeshData::getInstance().getData <vector<int> > ("M02M1");
}

//--------------------------------------------------------------------------//

void Interp::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  XYSh =
   PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

