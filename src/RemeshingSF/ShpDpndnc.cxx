// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "ShpDpndnc.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

ShpDpndnc::ShpDpndnc (double X,double Y,double USH,
                      double VSH,double XI,double YI,
                      double UI,double VI,double AI)
{
 x=X; y=Y; ush=USH; vsh=VSH;
 xi=XI; yi=YI; ui=UI; vi=VI; ai=AI;
}

//--------------------------------------------------------------------------//

ShpDpndnc::~ShpDpndnc()
{
}

//--------------------------------------------------------------------------//

unsigned ShpDpndnc::getDependence()
{
  callShpDpndnc();
  return shp_dependence;
}

//--------------------------------------------------------------------------//

void ShpDpndnc::callShpDpndnc ()
{
  double rl, umod, dt;
  double xx, yy, xxi, yyi;
  double dsh1, dsh2;

  //compute dt
  rl = sqrt(pow((x-xi),2)+pow((y-yi),2));
  umod = sqrt(ui*ui+vi*vi);
  dt = 0.01*rl/(ai+umod);
  //compute shock position after dt
  xx = x+ush*dt;
  yy = y+vsh*dt;
  // compute position of i-point after dt
  xxi = xi+ui*dt;
  yyi = yi+vi*dt;
  // compute distance between i-ipoin and
  // shock point
  dsh1 = sqrt(pow((x-xi),2)+pow((y-yi),2));
  // compute distance between i-ipoin and
  // shock point after dt=1
  dsh2 = pow((xx-xxi),2)+pow((yy-yyi),2);
  dsh2 = sqrt(dsh2)-ai*dt;

  // compute distances
  shp_dependence = 0;
  if(dsh2<dsh1) { shp_dependence = 1;}
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
