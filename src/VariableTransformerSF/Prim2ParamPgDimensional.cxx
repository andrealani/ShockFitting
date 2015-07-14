// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "VariableTransformerSF/Prim2ParamPgDimensional.hh"
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
   kinetic = (*zroe)(1,IPOIN) * (*zroe)(1,IPOIN) +
             (*zroe)(2,IPOIN) * (*zroe)(2,IPOIN);
   kinetic = kinetic*0.5;
   help = ReferenceInfo::getgam()/(ReferenceInfo::getgam()-1);
   h = help * ReferenceInfo::getRgas() * (*zroe)(3,IPOIN) + kinetic;

   (*zroe)(3,IPOIN) = sqrtr * (*zroe)(2,IPOIN) / ReferenceInfo::geturef();
   (*zroe)(2,IPOIN) = sqrtr * (*zroe)(1,IPOIN) / ReferenceInfo::geturef();
   (*zroe)(1,IPOIN) = sqrtr * h / (ReferenceInfo::geturef() *
                                   ReferenceInfo::geturef());
   (*zroe)(0,IPOIN) = sqrtr;

    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     (*XY)(I,IPOIN) = (*XY)(I,IPOIN)/ReferenceInfo::getLref(); }
  } 

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void Prim2ParamPgDimensional::transform(vector <double>& m_prim, 
					vector <double>& m_XY,
                                        vector <double>& m_zroe)
{
  double sqrtr;

  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);

  rho = m_prim.at(0) / m_prim.at(3) / ReferenceInfo::getRgas();
  sqrtr = sqrt(rho / ReferenceInfo::getrhoref());
  kinetic = m_prim.at(1) * m_prim.at(1) + m_prim.at(2) * m_prim.at(2);
  kinetic = kinetic*0.5;
  help = ReferenceInfo::getgam()/(ReferenceInfo::getgam()-1);
  h = help * ReferenceInfo::getRgas() * m_prim.at(3) + kinetic;

  m_zroe.at(3) = sqrtr * m_prim.at(2) / ReferenceInfo::geturef();
  m_zroe.at(2) = sqrtr * m_prim.at(1) / ReferenceInfo::geturef();
  m_zroe.at(1) = sqrtr * h / (ReferenceInfo::geturef() *
                                 ReferenceInfo::geturef());
  m_zroe.at(0) = sqrtr; 

  m_XY.at(0) = m_XY.at(0)/ReferenceInfo::getLref();
  m_XY.at(1) = m_XY.at(1)/ReferenceInfo::getLref();

}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
