// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Param2PrimPgDimensional.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
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

  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    rho = (*zroe)(0,IPOIN) * (*zroe)(0,IPOIN);
    kinetic = (*zroe)(2,IPOIN) * (*zroe)(2,IPOIN) + 
              (*zroe)(3,IPOIN) * (*zroe)(3,IPOIN);
    kinetic = kinetic * 0.5;
    h = (*zroe)(1,IPOIN)/(*zroe)(0,IPOIN);
    help = (ReferenceInfo::getgam()-1)/ReferenceInfo::getgam();
    pres = help * (rho * h - kinetic);
    u = (*zroe)(2,IPOIN)/(*zroe)(0,IPOIN);
    v = (*zroe)(3,IPOIN)/(*zroe)(0,IPOIN);

    (*zroe)(0,IPOIN) = pres * ReferenceInfo::getrhoref() *
                       ReferenceInfo::geturef() * ReferenceInfo::geturef();
    (*zroe)(1,IPOIN) = u * ReferenceInfo::geturef();
    (*zroe)(2,IPOIN) = v * ReferenceInfo::geturef();
    (*zroe)(3,IPOIN) = pres/rho *
                       (ReferenceInfo::geturef() * ReferenceInfo::geturef() /
                        ReferenceInfo::getRgas());

    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     (*XY)(I,IPOIN) = (*XY)(I,IPOIN)*ReferenceInfo::getLref(); }
  }

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void Param2PrimPgDimensional::transform(vector <double>* m_zroe,
                                        vector <double>* m_XY,
                                        vector <double>* m_prim)
{
  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);

  rho = m_zroe->at(0) * m_zroe->at(0);
  kinetic = m_zroe->at(2) * m_zroe->at(2) +
            m_zroe->at(3) * m_zroe->at(3);
  kinetic = kinetic * 0.5;
  h = m_zroe->at(1)/m_zroe->at(0);
  help = (ReferenceInfo::getgam()-1)/ReferenceInfo::getgam();
  pres = help * (rho * h - kinetic);
  u = m_zroe->at(2)/m_zroe->at(0);
  v = m_zroe->at(3)/m_zroe->at(0);

  m_prim->at(0) = pres * ReferenceInfo::getrhoref() *
                       ReferenceInfo::geturef() * ReferenceInfo::geturef();
  m_prim->at(1) = u * ReferenceInfo::geturef();
  m_prim->at(2) = v * ReferenceInfo::geturef();
  m_prim->at(3) = pres/rho *
                       (ReferenceInfo::geturef() * ReferenceInfo::geturef() /
                        ReferenceInfo::getRgas());

  for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
   m_XY->at(I) = m_XY->at(I)/ReferenceInfo::getLref();
  }
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
