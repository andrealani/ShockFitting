// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "MathTools/Ishel.hh"
#include "MathTools/Solg.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

Ishel::Ishel ()
{
}

//--------------------------------------------------------------------------//

Ishel::Ishel (std::vector <double> XC, std::vector <double> YC,
               double XS1, double XS2, double YS1, double YS2)
{
 xc = XC; yc = YC; xs1 = XS1; xs2 = XS2; ys1 = YS1; ys2 = YS2;
}

//--------------------------------------------------------------------------//

Ishel::~Ishel()
{
}

//--------------------------------------------------------------------------//

int Ishel::Ishel1 ()
{
  double eval;
  int ishel1=0;
  for (unsigned i=0; i<xc.size(); i++) {
   eval = (xs1-xc.at(i))*(ys1-ys2)-(ys1-yc.at(i))*(xs1-xs2);
   ishel1 = ishel1 +copysign(1,eval);}
  ishel1 = abs(ishel1)/3;
  return ishel1;
}

//--------------------------------------------------------------------------//

int Ishel::Ishel2 ()
{
  Array2D <double> a;
  vector <double> b;
  vector <double> x;
  double xi, yi, rlsh2, rl2;
  a.resize(2,2);
  b.resize(2);
  x.resize(2);

  // write the equation of the straight line passing for the two shock points
  a(0,0) = ys2-ys1;
  a(0,1) = xs1-xs2;
  b.at(0) = xs1*(ys2-ys1)+ys1*(xs1-xs2);

  // write the equation of the straight line perpendicular to previous shock line
  // and passing for the vertex node
  a(1,0)=(yc.at(1)-yc.at(0));
  a(1,1)=(xc.at(0)-xc.at(1));
  b.at(1)=xc.at(0)*(yc.at(1)-yc.at(0))+yc.at(0)*(xc.at(0)-xc.at(1));

  //solve the linear system and find the intersection point coordinates
  Solg <double> S;
  x = S.callSolg(a,b);
  xi = x.at(0);
  yi = x.at(1);

  //compute the distance between the shock points and the sum of  distances
  // of the  intersection point from the both shock points
  rlsh2=pow((xs1-xs2),2)+pow((ys1-ys2),2);
  rl2 = pow((xs1-xi),2)+pow((ys1-yi ),2)+pow((xs2-xi ),2)+pow((ys2-yi ),2);

  //if the distance between the shock points is equal to the sum of distances
  //of the  intersection point from the both shock points the intersection
  //point is enclosed between the two shock point
  if (rlsh2 >= rl2) return 0;

  //repeat the same algorithm for other triangle vertices
  a(0,0)=(ys2-ys1);
  a(0,1)=(xs1-xs2);
  b.at(0)=xs1*(ys2-ys1)+ys1*(xs1-xs2);

  a(1,2)=(yc.at(2)-yc.at(0));
  a(1,1)=(xc.at(0)-xc.at(2));
  b.at(1)=xc.at(0)*(yc.at(2)-yc.at(0))+yc.at(0)*(xc.at(0)-xc.at(2));

  x = S.callSolg(a,b);
  xi = x.at(0);
  yi = x.at(1);

  rlsh2=pow((xs1-xs2),2)+pow((ys1-ys2),2);
  rl2  =pow((xs1-xi),2)+pow((ys1-yi ),2)+pow((xs2-xi ),2)+pow((ys2-yi ),2);
  if (rlsh2 >=rl2) return 0;

  //
  a(0,0)=(ys2-ys1);
  a(0,1)=(xs1-xs2);
  b.at(0)=xs1*(ys2-ys1)+ys1*(xs1-xs2);

  a(1,0)=(yc.at(2)-yc.at(1));
  a(1,1)=(xc.at(1)-xc.at(2));
  b.at(1)=xc.at(1)*(yc.at(2)-yc.at(1))+yc.at(1)*(xc.at(1)-xc.at(2));

  x = S.callSolg(a,b);

  xi=x.at(0);
  yi=x.at(1);

  rlsh2=pow((xs1-xs2),2)+pow((ys1-ys2),2);
  rl2  =pow((xs1-xi),2)+pow((ys1-yi),2)+pow((xs2-xi),2)+pow((ys2-yi),2);
  if(rlsh2 >= rl2) return 0;

  return 1;
}

