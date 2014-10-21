// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/FindBEdg.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"

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
  setPhysicsData();

  setAddress();
  unsigned ibfac=0;

  while (ibfac < *nbfac) {
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

     if (0<=s && s<=1) { return ibfac;} 
    } // if abs(tline)<=tollerance
   ibfac++;
   } // if bndfac(2,ibfac)>0

   else {ibfac++;}
  } // while

  return -1;
}

//--------------------------------------------------------------------------//

void FindBEdg::setAddress()
{
  unsigned start;
  start = 0; 
  XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start)); 
}

//--------------------------------------------------------------------------//

void FindBEdg::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  coor = MeshData::getInstance().getData <std::vector<double> > ("COOR");
  bndfac = MeshData::getInstance().getData <Array2D<int> > ("BNDFAC"); 
}

//--------------------------------------------------------------------------//

void FindBEdg::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
