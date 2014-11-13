// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/CoDc.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/PhysicsData.hh"
#include "MathTools/Solg.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

CoDc::CoDc()
{
}

//----------------------------------------------------------------------------//

CoDc::~CoDc()
{
}

//----------------------------------------------------------------------------//

void CoDc::callCoDc(vector<double> xd, vector<double> xu)
{
  setPhysicsData();

  gm1 = (*gam)-1;

  x1.resize(xd.size()); x2.resize(xu.size());
  x1 = xd; x2 = xu; 
  W=0;

  // resize local arrays
  yn1.resize(7);
  yn.resize(7);
  b.resize(7);
  g.resize(7,7);
  dyn.resize(7);

  // compute riemann invariants
  R1 = sqrt( (*gam) * x1.at(1)/x1.at(0) ) + gm1/2 * x1.at(2);
  R2 = sqrt( (*gam) * x2.at(1)/x2.at(0) ) - gm1/2 * x2.at(2);
  S1 = pow( (x1.at(1)/x1.at(0)) , (*gam) );
  S2 = pow ( (x2.at(1)/x2.at(0)) , (*gam) );

  // initialize unknows vector
  yn1.at(0) = x1.at(0);
  yn1.at(1) = x1.at(1);
  yn1.at(2) = x1.at(2);
  yn1.at(3) = x2.at(0);
  yn1.at(4) = x2.at(1);
  yn1.at(5) = x2.at(2);
  yn1.at(6) = W;

  dum = 1;

  while (dum > pow(10,(-7))){
   for (unsigned i=0; i<7; i++) {yn.at(i) = yn1.at(i);}

   // compute the Jacobian matrix
   for (unsigned i=0; i<7; i++) {
    b.at(i) = fDc(i,yn);
    for (unsigned j=0; j<7; j++) {
     for (unsigned k=0; k<7; k++) { yn1.at(k) = yn.at(k); }
     dum1 = fDc(i,yn);
     dyn1 = abs(yn.at(j))*0.01;
     if (dyn1 <= pow(10,(-7))) {dyn1 = pow(10,(-7));}
     yn1.at(j) = yn.at(j)-dyn1;
     dum1 = fDc(i,yn1);
     yn1.at(j) = yn.at(j)+dyn1;
     dum2 = fDc(i,yn1);
     g(i,j)=(dum2-dum1)/(2.0*dyn1);
    }
   }

    Solg <double> System;
    dyn = System.callSolg(g,b);
    for (unsigned i=0; i<7; i++) { yn1.at(i) = yn.at(i)-dyn.at(i); } 

    dum =0;
    for (unsigned i=0; i<7; i++) { dum = dum+abs(yn1.at(i)-yn.at(i)); }
  }

  W = yn1.at(6);
  x1.at(0) = yn1.at(0);
  x1.at(1) = yn1.at(1);
  x1.at(2) = yn1.at(2);
  x2.at(0) = yn1.at(0);
  x2.at(1) = yn1.at(1);
  x2.at(2) = yn1.at(2);
}

//----------------------------------------------------------------------------//

double CoDc::fDc(unsigned index, vector<double> y)
{
  if (index==0) {
   return sqrt( (*gam) * y.at(1)/y.at(0)) + gm1/2 * y.at(2) - R1;}
  else if (index==1) {
   return pow( (y.at(1)/y.at(0)), (*gam) ) - S1;}
  else if (index==2) {
   return sqrt( (*gam) * y.at(4)/y.at(3)) - gm1/2 * y.at(5) - R2;}
  else if (index==3) {
   return pow( (y.at(4)/y.at(3)) , (*gam)) - S2;}
  else if (index==4) {
   return y.at(1) - y.at(4);}
  else if (index==5) {
   return y.at(2) - y.at(5);}
  else if (index==6) {
   return y.at(6) - y.at(2);}
  else {
   cout << "ComputeStateDps::error => index error in private member fDc ";
   cout << "of CoDc\n";
   exit(1);}

}

//----------------------------------------------------------------------------//

void CoDc::setPhysicsData()
{
  if(ChemicalInfo::getModel()=="PG") {
   *gam = PhysicsInfo::getGam();
  }
  else if (ChemicalInfo::getModel() == "Cneq" ||
           ChemicalInfo::getModel() == "TCneq") {
   gam = PhysicsData::getInstance().getData <double >("GREF");
  }
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting


