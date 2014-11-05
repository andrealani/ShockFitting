// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/CoUTP.hh"
#include "Framework/PhysicsData.hh"
#include "MathTools/Solg.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

CoUTP::CoUTP()
{
}

//----------------------------------------------------------------------------//

CoUTP::~CoUTP()
{
}

//----------------------------------------------------------------------------//

void CoUTP::callCoUTP(vector<double> y, double r14, double r23,
                      double dxR14, double dyR14, double WWS)
{
  R14 = r14; R23 = r23; 
  DXR14 = dxR14; DYR14 = dyR14; 
  unSh = WWS;

  const unsigned nn=20;
  unsigned icont;
  double sx14, taux14, tauy14, nx14, ny14;
  vector <double> bb(nn);
  Array2D <double> g(nn,nn);

  yn.resize(nn);
  a.resize(11,nn);
  b.resize(nn);

  ifail = false;

  logfile.Open(getClassName());

  // set Array2D a and vector b
  for(unsigned i=0; i<nn-(2*4+1); i++) {
   for(unsigned j=0; j<nn; j++) { a(i,j) = 0.0; }
   b.at(i) = 0.0;
  }
  // zone (1): known zone between incident shock and mach stem
  a(0,0) = 1.0;
  a(1,1) = 1.0;
  a(2,2) = 1.0;
  a(3,3) = 1.0;
  b.at(0) = y.at(0);
  b.at(1) = y.at(1);
  b.at(2) = y.at(2);
  b.at(3) = y.at(3);
  // zone (2): known zone between incident shock and reflected shock
  a(4,4) = 1.0;
  a(5,5) = 1.0;
  a(6,6) = 1.0;
  a(7,7) = 1.0;
  b.at(4) = y.at(4);
  b.at(5) = y.at(5);
  b.at(6) = y.at(6);
  b.at(7) = y.at(7);
  // known slope of the incident shock
  a(8,16) = 1.0;
  b.at(8) = y.at(16);

  // initialize unknowns vector
  for(unsigned i=0; i<nn; i++) { yn1.at(i) = y.at(i); }

  icont = 0;
ten:

  for(unsigned i=0; i<nn; i++) {
   yn.at(i) = yn1.at(i);
   bb.at(i) = fUTP(yn1,i);
  }
  for(unsigned i=0; i<nn; i++) {
   for(unsigned j=0; j<nn; j++) {
    for(unsigned k=0; k<nn; k++) { yn1.at(k)=yn.at(k); }
    dyn1 = abs(yn1.at(j))*.001;
    if(dyn1<pow(10,-7)) { dyn1 = pow(10,-7); }
    yn1.at(j) = yn.at(j)+dyn1;
    dum2 = fUTP(yn1,i);
    yn1.at(j) = yn.at(j)-dyn1;
    dum1 = fUTP(yn1,i);
    g(i,j)=(dum2-dum1)/(2*dyn1);
   }
  }

  Solg <double> System;
  dyn = System.callSolg(g,b);

  for(unsigned i=0; i<nn; i++) { yn1.at(i)=yn.at(i) - 0.5 * dyn.at(i); }

  dum = 0;
  for(unsigned i=0; i<nn; i++) { dum = dum + abs(yn1.at(i)-yn.at(i)); }

  ++icont;
  logfile ("Conv ---> ", dum, " ", icont, "\n");

  // compute R14 with the new value of n14;
  sx14 = yn1.at(18);
  taux14 = cos(sx14);
  tauy14 = sin(sx14);
  nx14 = -tauy14;
  ny14 = taux14;

  if (icont==1) { dumold = dum; goto ten; }

  if (dum>dumold && icont>1) { ifail = true; return; }

  if (dum>pow(10,-7)) { dumold = dum; goto ten; }

  logfile("Initial and final state:\n");
  for(unsigned i=0; i<nn; i++) {
   logfile(i," ",y.at(i)," ", yn.at(i), "\n");
  }

  logfile.Close();
}

//----------------------------------------------------------------------------//

double CoUTP::fUTP(std::vector<double> y, unsigned index)
{
  double ro1, p1, u1, v1;
  double ro2, p2, u2, v2;
  double sx12, sx14, sx23;
  double wsh, wn, wt;
  double taux12, tauy12, nx12, ny12;
  double taux14, tauy14, nx14, ny14;
  double taux23, tauy23, nx23, ny23;
  double un1, un2, ut1, ut2;


  if(index >= 0 && index < 4) {
   ro1 = y.at(4); p1 = y.at(5); u1 = y.at(6); v1 = y.at(7);
   ro2 = y.at(8); p2 = y.at(9); u2 = y.at(10); v2 = y.at(11);
   sx12 = y.at(16); sx23 = y.at(17);
   wsh = y.at(19); 
   // compute unit vector tangential to the incident shock
   taux12 = cos(sx12);
   tauy12 = sin(sx12);
   nx12 = -tauy12;
   ny12 = taux12; 
   // compute unit vector tangential to the reflected shock
   taux23 = cos(sx23);
   tauy23 = sin(sx23);
   nx23 = -tauy23;
   ny23 = taux23;
   // compute normal and tangential components of the speed
   un1 = u1*nx23+v1*ny23;
   un2 = u2*nx23+v2*ny23;
   ut1 = u1*taux23+v1*tauy23;
   ut2 = u2*taux23+v2*tauy23;
   // compute normal speed of the triple point
   wn = wsh * (taux12*nx23+tauy12*ny23);
   // compute tangential speed of the triple point
   wt = wsh * (taux12*taux23+tauy12*tauy23) + 
        unSh * (nx12  *taux23+ny12  *tauy23);
   if      (index==0)  {
    return ro1*un1-ro2*un2-wn*(ro1-ro2); }
   else if (index==1)  {
    return p1+ro1*pow((un1-wn),2)-p2-ro2*pow((un2-wn),2); }
   else if (index==2)  {
    return ut1-ut2; }
   else if (index==3)  {
    return (*gam)/((*gam)-1.0) * p1/ro1 + 0.5 * pow((un1-wn),2) -
           (*gam)/((*gam)-1.0) * p2/ro2 - 0.5 * pow((un2-wn),2); }
   else                {
    cout << "CoUTP::error => wrong index, here it must be < 4\n";
    exit(1);  }
  }

  else if (index >= 4 && index < 9) {
   ro1 = y.at(0); p1 = y.at(1); u1 = y.at(2); v1 = y.at(3);
   ro2 = y.at(12); p2 = y.at(13); u2 = y.at(14); v2 = y.at(15);
   sx12 = y.at(16); sx14 = y.at(18);
   wsh = y.at(19);
   // compute unit vector tangential to the incident shock
   taux12 = cos(sx12);
   tauy12 = sin(sx12);
   nx12 = -tauy12;
   ny12 = taux12;
   // compute unit vector tangential to the reflected shock
   taux14 = cos(sx14);
   tauy14 = sin(sx14);
   nx14 = -tauy14;
   ny14 = taux14;
   // compute normal and tangential components of the speed
   un1 = u1*nx14+v1*ny14;
   un2 = u2*nx14+v2*ny14;
   ut1 = u1*taux14+v1*tauy14;
   ut2 = u2*taux14+v2*tauy14;
   // compute normal speed of the triple point
   wn = wsh * (taux12*nx14+tauy12*ny14);
   // compute tangential speed of the triple point
   wt = wsh * (taux12*taux14+tauy12*tauy14)+
        unSh * (nx12  *taux14+ny12  *tauy14);
   if      (index==4) {
    return ro1*un1-ro2*un2-wn*(ro1-ro2); }
   else if (index==5) {
    return p1+ro1*pow((un1-wn),2)-p2-ro2*pow((un2-wn),2); }
   else if (index==6) {
    return ut1-ut2; }
   else if (index==7) {
    return (*gam)/((*gam)-1.0) * p1/ro1 + 0.5 * pow((un1-wn),2) -
           (*gam)/((*gam)-1.0) * p2/ro2 - 0.5 * pow((un2-wn),2); }
   else if (index==8) {
    ro2 = y.at(12); p2 = y.at(13); u2 = y.at(14); v2 = y.at(15);
    sx14 = y.at(18);
    // compute unit vector tangential to the reflected shock
    taux14 = cos(sx14);
    tauy14 = sin(sx14);
    nx14 = -tauy14;
    ny14 = taux14;
    un2 = u2*DXR14+v2*DYR14;
    return sqrt((*gam) * p2/ro2) + delta * un2 - R14;
   }
   else                { 
    cout << "CoUTP::error => wrong index, here it must be included between 5 and 9\n";
    exit(1);  }
  }

  else if (index==9 || index==10) {
   ro1 = y.at(8); p1 = y.at(9); u1 = y.at(10); v1 = y.at(11);
   ro2 = y.at(12); p2 = y.at(13); u2 = y.at(14); v2 = y.at(15);
   if      (index==9)  { return p1-p2; }
   else if (index==10) { return (u1*v2-v1*u2); }
   else { 
    cout << "CoUTP::error => wrong index, here it must be equal to 9 or 10\n";
    exit(1); }
  }

  else if(index >= 11 && index<20) {
   double futp=0;
   unsigned ii = index-(2*4+1+1+1);
   for(unsigned j=0; j<20; j++) { futp=futp+a(ii,j)*y.at(j); }
   return futp-b.at(ii);
  }
  else  {
   cout << "CoUTP::error => wrong index, it must be included between 0 and 19\n";
    exit(1);  }

} 

//----------------------------------------------------------------------------//

void CoUTP::setPhysicsData()
{
  gam = PhysicsData::getInstance().getData <double> ("GAM");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
