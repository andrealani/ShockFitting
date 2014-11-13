// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Prim2ParamPgDimensional.hh"
#include "Framework/Log.hh"
#include "Framework/ReferenceInfo.hh"
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

  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
   rho = (*zroe)(0,IPOIN) / (*zroe)(3,IPOIN) / ReferenceInfo::getRgas();
   sqrtr = sqrt(rho / ReferenceInfo::getrhoref());
   kinetic = pow((*zroe)(1,IPOIN),2) + pow((*zroe)(2,IPOIN),2);
   kinetic = kinetic*0.5;
   help = ReferenceInfo::getgam()/(ReferenceInfo::getgam()-1);
   h = help * ReferenceInfo::getRgas() * (*zroe)(3,IPOIN) + kinetic;

   (*zroe)(3,IPOIN) = sqrtr * (*zroe)(2,IPOIN) / ReferenceInfo::geturef();
   (*zroe)(2,IPOIN) = sqrtr * (*zroe)(1,IPOIN) / ReferenceInfo::geturef();
   (*zroe)(1,IPOIN) = sqrtr * h / pow(ReferenceInfo::geturef(),2);
   (*zroe)(0,IPOIN) = sqrtr;
  } 
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
