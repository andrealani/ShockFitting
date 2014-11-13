// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/CoShock.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/PhysicsData.hh"
#include "MathTools/Solg.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

CoShock::CoShock()
{
}

//----------------------------------------------------------------------------//

CoShock::~CoShock()
{
}

//----------------------------------------------------------------------------//

void CoShock::callCoShock(vector<double> xd, vector<double> xu, double R )
{
  setPhysicsData();

  x1.resize(xd.size()); x2.resize(xu.size());
  x1 = xd; x2 = xu; R2 = R;
  W=0;

  // resize local arrays
  yn1.resize(4);
  yn.resize(4);
  b.resize(4);
  g.resize(4,4);
  dyn.resize(4);

  // initialize unknows vectors
  yn1.at(0) = x1.at(0);
  yn1.at(1) = x1.at(1);
  yn1.at(2) = x1.at(2);
  yn1.at(3) = W;

  dum = 1;

  while (dum > pow(10,(-7))) {
   for (unsigned i=0; i<4; i++) { yn.at(i)=yn1.at(i); }

   // compute the jacobian matrix
   for (unsigned i=0; i<4; i++) {
    b.at(i) = fShock(i, yn);
    for (unsigned j=0; j<4; j++) {
     for (unsigned k=0; k<4; k++) { yn1.at(k) = yn.at(k); }
     dyn1 = abs(yn1.at(j))*0.01;
     if (dyn1 <= pow(10,(-8))) {dyn1 = pow(10,(-7));}
     yn1.at(j) = yn.at(j)-dyn1;
     dum1 = fShock(i, yn1);
     yn1.at(j) = yn.at(j)+dyn1;
     dum2 = fShock(i, yn1);
     g(i,j) = (dum2-dum1)/(2*dyn1);
    }
   }

   Solg <double> System;
   dyn = System.callSolg(g,b);

   for(unsigned i=0; i<4; i++) { yn1.at(i) = yn.at(i) -0.2*dyn.at(i); }

   dum=0;
   for (unsigned i=0; i<4; i++) { dum = dum +abs(yn1.at(i)-yn.at(i)); }
  }

  W = yn1.at(3);
  x1.at(0) = yn1.at(0);
  x1.at(1) = yn1.at(1);
  x1.at(2) = yn1.at(2);
}

//----------------------------------------------------------------------------//

double CoShock::fShock (unsigned index, vector <double> y)
{
  double  gm1 = (*gam)-1;
  if (index==0){
   return y.at(0)*(y.at(2) - y.at(3)) - x2.at(0) * (x2.at(2) - y.at(3));}
  else if (index==1) {
   return y.at(1) + y.at(0) * pow((y.at(2) - y.at(3)),2)-
          x2.at(1) - x2.at(0) * pow((x2.at(2) - y.at(3)),2);}
  else if (index==2) {
   return (*gam)/gm1 * y.at(1)/y.at(0) + 0.5 * pow((y.at(2) - y.at(3)),2)
          - (*gam)/gm1 * x2.at(1)/x2.at(0) - 0.5 * pow((x2.at(2) - y.at(3)),2);} 
  else if (index==3) {
   return sqrt((*gam) * y.at(1)/y.at(0)) + gm1/2 * y.at(2) - R2;}
  else {
   cout << "ComputeStateDps::error => index error in private member fShock ";
   cout << "of class CoShock\n";
   exit(1);}
}

//----------------------------------------------------------------------------//

void CoShock::setPhysicsData()
{
  if(ChemicalInfo::getModel() == "PG") {
   gam = PhysicsData::getInstance().getData <double >("GAM");
  }
  else if (ChemicalInfo::getModel() == "Cneq" ||
           ChemicalInfo::getModel() =="TCneq") {
   gam = PhysicsData::getInstance().getData <double >("GREF");
  }
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
