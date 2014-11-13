// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include <vector>
#include "Framework/MeshData.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Solg.hh"
#include "RemeshingSF/Rdshp.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Rdshp::Rdshp ()
{
}

//--------------------------------------------------------------------------//

Rdshp::~Rdshp()
{
}

//--------------------------------------------------------------------------//

double Rdshp::getRdshp(double XC, double YC, double XS1,
                       double YS1, double XS2, double YS2) 
{
  resetValues(XC,YC,XS1,YS1,XS2,YS2);
  callRdshp();
  return rdshp;
}

//--------------------------------------------------------------------------//

void Rdshp::callRdshp ()

{
  Array2D <double> a (2,2);
  vector <double> b(2);
  vector <double> x(2);
  double xi, yi;
  double rlsh2, rlsh3, rl2;

  rdshp = -1;
  a(0,0) = (ys2-ys1);
  a(0,1) = (xs1-xs2);
  b.at(0) = xs1*(ys2-ys1)+ys1*(xs1-xs2);
  a(1,0) = a(0,1);
  a(1,1) = -a(0,0);
  b.at(1) = a(1,0)*xc+a(1,1)*yc;

  Solg <double> System;
  x = System.callSolg(a,b);
  xi = x.at(0);
  yi = x.at(1);

  rlsh2 = (xs1-xs2)*(xs1-xs2)+(ys1-ys2)*(ys1-ys2);
  rlsh3=((1.0+(MeshData::getInstance().getSNDMIN())) * 
         (1.0+MeshData::getInstance().getSNDMIN()) +
         (MeshData::getInstance().getSNDMIN()) * 
         (MeshData::getInstance().getSNDMIN()))*rlsh2;
  rl2 = (xs1-xi)*(xs1-xi)+(ys1-yi)*(ys1-yi)+
        (xs2-xi)*(xs2-xi)+(ys2-yi)*(ys2-yi);
  if (rlsh3 < rl2) {return;}
  rdshp = (xi-xc)*(xi-xc)+(yi-yc)*(yi-yc);
  rdshp = rdshp/rlsh2;
  rdshp = sqrt(rdshp);
  return;
}

//--------------------------------------------------------------------------//

void Rdshp::resetValues (double XC, double YC, double XS1,
                         double YS1, double XS2, double YS2)
{
  xc=XC, yc=YC; xs1=XS1; ys1=YS1; xs2=XS2; ys2=YS2; 
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

