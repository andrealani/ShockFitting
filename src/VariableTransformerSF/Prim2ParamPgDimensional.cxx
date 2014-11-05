// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Prim2ParamPgDimensional.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Prim2ParamPgDimensional, VariableTransformer>
 prim2ParamPgDimensionalProv("Prim2ParamPgDimensional");

//--------------------------------------------------------------------------//

Prim2ParamPgDimensional::Prim2ParamPgDimensional(const std::string& objectName) :
 Prim2Param(objectName)
{
}

//--------------------------------------------------------------------------//

Prim2ParamPgDimensional::~Prim2ParamPgDimensional()
{
}

//--------------------------------------------------------------------------//

void Prim2ParamPgDimensional::setup()
{
  LogToScreen(VERBOSE, "Prim2ParamPgDimensional::setup() => start \n");

  LogToScreen(VERBOSE, "Prim2ParamPgDimensional::setup() => end \n");
}

//--------------------------------------------------------------------------//

void Prim2ParamPgDimensional::unsetup()
{
  LogToScreen(VERBOSE, "Prim2ParamPgDimensional::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Prim2ParamPgDimensional::transform()
{
  LogToScreen(INFO, "Prim2ParamPgDimensional::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  double sqrtr;

  (*rhoref) = (*pref) / (*Tref) / (*Rgas);

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   rho = (*zroe)(0,IPOIN) / (*zroe)(3,IPOIN) / (*Rgas);
   sqrtr = sqrt(rho / (*rhoref));
   kinetic = pow((*zroe)(1,IPOIN),2) + pow((*zroe)(2,IPOIN),2);
   kinetic = kinetic*0.5;
   help = (*gam)/((*gam)-1);
   h = help * (*Rgas) * (*zroe)(3,IPOIN) + kinetic;

   (*zroe)(3,IPOIN) = sqrtr * (*zroe)(2,IPOIN) / (*uref);
   (*zroe)(2,IPOIN) = sqrtr * (*zroe)(1,IPOIN) / (*uref);
   (*zroe)(1,IPOIN) = sqrtr * h / pow((*uref),2);
   (*zroe)(0,IPOIN) = sqrtr;
  } 
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
