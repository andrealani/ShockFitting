// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/ComputeTv.hh"
#include "Framework/ChemicalConsts.hh"
#include "Framework/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

ComputeTv::ComputeTv()
{
  F=0; Fd=0;
  tollerance = pow(10,-12); tolleranceD = pow(10,-15);
  KMAX = 100;
}

//--------------------------------------------------------------------------//

ComputeTv::~ComputeTv()
{
}

//--------------------------------------------------------------------------//

void ComputeTv::callComputeTv(double Ev, vector<double> ALPHA, 
                              vector<double> Temp)
{
  T = Temp; ev = Ev; alpha = ALPHA;

  double dT, delta, Td;
  unsigned K;

  setPhysicsData();

  T.at(1)=T.at(0);
  dT=0.001*T.at(0);

  funct();

  delta = abs(2*dT/T.at(1));

  K = 0;

  while(abs(delta)>tolleranceD && K<=KMAX) {
   funct();

   Td = T.at(1)-F/Fd;
   dT = Td-T.at(1);
   T.at(1) = Td;
   delta = abs(dT/T.at(1));
   ++K;
  }
}

//--------------------------------------------------------------------------//

void ComputeTv::funct()
{
  double thelp;
  std::vector<double> evs(*nsp);
  std::vector<double> cvs(*nsp);

  F =0;
  Fd = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
    thelp = thev->at(ISP)/T.at(1);

    if(typemol->at(ISP)=="A")       { evs.at(ISP) = 0;
                                      cvs.at(ISP) = 0; }
    else if (typemol->at(ISP)=="B") {
     evs.at(ISP) = ChemicalConsts::Rgp()/mm->at(ISP) *
                   thev->at(ISP)/(exp(thelp)-1);
     cvs.at(ISP) = ChemicalConsts::Rgp()/mm->at(ISP) *
                   pow(thelp,2)*exp(thelp)/pow((exp(thelp)-1),2);
    }
    else                             { evs.at(ISP) = 0;
                                       cvs.at(ISP) = 0; } 

    F = F + alpha.at(ISP)*evs.at(ISP);
    Fd = Fd + alpha.at(ISP)*cvs.at(ISP);
  }
  F = ev - F;
  Fd = (-1)*Fd;
}

//--------------------------------------------------------------------------//

void ComputeTv::setPhysicsData()
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  thev = PhysicsData::getInstance().getData <std::vector<double> > ("THEV");
  mm = PhysicsData::getInstance().getData <std::vector<double> > ("MM");
  typemol =
    PhysicsData::getInstance().getData <std::vector<std::string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
