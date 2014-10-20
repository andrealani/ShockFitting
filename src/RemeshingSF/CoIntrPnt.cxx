// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/CoIntrPnt.hh"
#include "MathTools/Solg.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

CoIntrPnt::CoIntrPnt()
{
}

//--------------------------------------------------------------------------//

CoIntrPnt::~CoIntrPnt()
{
}

//--------------------------------------------------------------------------//

double CoIntrPnt::getX() const { return x.at(0);}

//--------------------------------------------------------------------------//

double CoIntrPnt::getY() const { return x.at(1);}

//--------------------------------------------------------------------------//

void CoIntrPnt::callCoIntrPnt(vector <double> xc, vector <double> yc,
                              vector <double> xs, vector <double> ys)
{
  x.resize(2);
  Array2D <double> a(2,2);
  vector <double> b(2);

  a(0,0) = (ys.at(1)-ys.at(0));
  a(0,1) = (xs.at(0)-xs.at(1));
  b.at(0) = xs.at(0) * (ys.at(1)-ys.at(0)) + ys.at(0) * (xs.at(0)-xs.at(1));

  a(1,0) = yc.at(1)-yc.at(0);
  a(1,1) = xc.at(0)-xc.at(1);
  b.at(1) = xc.at(0) * (yc.at(1)-yc.at(0)) + yc.at(0) * (xc.at(0)-xc.at(1));

  Solg <double> S;
  x = S.callSolg(a,b);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
