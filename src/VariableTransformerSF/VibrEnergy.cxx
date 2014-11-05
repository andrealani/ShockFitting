// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/VibrEnergy.hh"
#include "Framework/ChemicalConsts.hh"
#include "Framework/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

VibrEnergy::VibrEnergy()
{
}

//--------------------------------------------------------------------------//

VibrEnergy::~VibrEnergy()
{
}

//--------------------------------------------------------------------------//

void VibrEnergy::callVibrEnergy(double Tv, vector<double> alpha)
{
  vector<double> evs(*nsp);

  setPhysicsData();

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   if      (typemol->at(ISP)=="A") {
    evs.at(ISP) = 0;
   }
   else if (typemol->at(ISP)=="B") {
    evs.at(ISP) = ChemicalConsts::Rgp() / mm->at(ISP) * thev->at(ISP) /
                  (exp(thev->at(ISP)/Tv)-1);
   }
   else                            {
    cout << "VibrEnergy::error => not implemented for ";
    cout << typemol->at(ISP) << " typemol \n";
   }
  }

  ev = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   ev = ev + alpha.at(ISP) * evs.at(ISP);
  }
}

//--------------------------------------------------------------------------//

void VibrEnergy::setPhysicsData()
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  thev = PhysicsData::getInstance().getData <vector<double> > ("THEV");
  mm = PhysicsData::getInstance().getData <vector<double> > ("MM");
  typemol =
    PhysicsData::getInstance().getData <vector<string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
