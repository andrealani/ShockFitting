// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Param2PrimPgDimensional.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Param2PrimPgDimensional, VariableTransformer>
 param2PrimPgDimensionalProv("Param2PrimPgDimensional");

//--------------------------------------------------------------------------//

Param2PrimPgDimensional::Param2PrimPgDimensional(const std::string& objectName) :
 Param2Prim(objectName)
{
}

//--------------------------------------------------------------------------//

Param2PrimPgDimensional::~Param2PrimPgDimensional()
{
}

//--------------------------------------------------------------------------//

void Param2PrimPgDimensional::setup()
{
  LogToScreen(VERBOSE, "Param2PrimPgDimensional::setup() => start \n");

  LogToScreen(VERBOSE, "Param2PrimPgDimensional::setup() => end \n");
}

//--------------------------------------------------------------------------//

void Param2PrimPgDimensional::unsetup()
{
  LogToScreen(VERBOSE, "Param2PrimPgDimensional::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Param2PrimPgDimensional::transform()
{
  LogToScreen(INFO, "Param2PrimPgDimensional::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  (*rhoref) = (*pref) / (*Tref) / (*Rgas);

  for(unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
    rho = (*v_Zroe)(0,IPOIN) * (*v_Zroe)(0,IPOIN);
    kinetic = pow((*v_Zroe)(2,IPOIN),2)+pow((*v_Zroe)(3,IPOIN),2);
    kinetic = kinetic * 0.5;
    h = (*v_Zroe)(1,IPOIN)/(*v_Zroe)(0,IPOIN);
    help = ((*gam)-1)/(*gam);
    pres = help * (rho * h - kinetic);
    u = (*v_Zroe)(2,IPOIN)/(*v_Zroe)(0,IPOIN);
    v = (*v_Zroe)(3,IPOIN)/(*v_Zroe)(0,IPOIN);

    (*v_Zroe)(0,IPOIN) = pres * (*rhoref) * (*uref) * (*uref);
    (*v_Zroe)(1,IPOIN) = u * (*uref);
    (*v_Zroe)(2,IPOIN) = v * (*uref);
    (*v_Zroe)(3,IPOIN) = pres/rho *
                         ((*uref) * (*uref) / (*Rgas));
  }
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
