// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/FindBEdg.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

FindBEdg::FindBEdg()
{
}

//--------------------------------------------------------------------------//

FindBEdg::~FindBEdg()
{
}

//--------------------------------------------------------------------------//

int FindBEdg::getBEdg(double xsh, double ysh)
{
  double x1, y1, x2, y2, tline;
  const double tollerance = 1e-10;

  setMeshData();

  setAddress();

  unsigned ibfac=0;

  while (ibfac < nbfac->at(0)) {
   int ibndfac1 = (*bndfac)(0,ibfac);
   int ibndfac2 = (*bndfac)(1,ibfac);


   if((*bndfac)(2,ibfac)>0) {

    x1 = (*XY)(0,ibndfac1-1); // c++ indeces start from 0
    y1 = (*XY)(1,ibndfac1-1); // c++ indeces start from 0
    x2 = (*XY)(0,ibndfac2-1); // c++ indeces start from 0
    y2 = (*XY)(1,ibndfac2-1); // c++ indeces start from 0
    tline = (ysh-y1)*(x2-x1)-(xsh-x1)*(y2-y1);

    if(std::abs(tline) <= tollerance) {
     if (std::abs(y2-y1)>std::abs(x2-x1)) { s = (ysh-y1)/(y2-y1); }
     else                              	{ s = (xsh-x1)/(x2-x1); }

     if (0<=s && s<=1) {   delete bndfac; delete XY; return ibfac;} 
    } // if abs(tline)<=tollerance
   ibfac++;
   } // if bndfac(2,ibfac)>0

   else {ibfac++;}
  } // while

  // de-allocate dynamic arrays
  freeArray();

  return -1;
}

//--------------------------------------------------------------------------//

void FindBEdg::freeArray()
{
  delete bndfac; delete XY;
}

//--------------------------------------------------------------------------//

void FindBEdg::setAddress()
{
  unsigned start = 0;
  unsigned totsize = nbfac->at(0) + 2 * 
                     PhysicsInfo::getnbShMax() *
                     PhysicsInfo::getnbShEdgesMax();
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             npoin->at(0),
                             &coorVect->at(start)); 
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
}

//--------------------------------------------------------------------------//

void FindBEdg::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> > ("BNDFAC"); 
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
