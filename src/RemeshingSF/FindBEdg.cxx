// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

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
  for (unsigned ibfac=0; ibfac<(*nbfac); ibfac++) {
   while((*bndfac)(2,ibfac)>0) {
    int ibndfac1 = (*bndfac)(0,ibfac-1); //c++indeces start from 0
    int ibndfac2 = (*bndfac)(1,ibfac-1); //c++indeces start from 0
    x1 = (*XY)(0,ibndfac1);
    y1 = (*XY)(1,ibndfac1);
    x2 = (*XY)(0,ibndfac2);
    y2 = (*XY)(1,ibndfac2);
    tline = (ysh-y1)*(x2-x1)-(xsh-x1)*(y2-y1);
    if(abs(tline)<=tollerance) {
     if (abs(y2-y1)>abs(x2-x1)) { s = (ysh-y1)/(y2-y1); }
     else			{ s = (xsh-x1)/(x2-x1); }
     if (0<=s && s<=1) { return ibfac;} 
    } // if abs(tline)<=tollerance
   } // while bndfac(2,ibfac)>0
  } // for
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
