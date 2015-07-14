// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "VariableTransformerSF/Prim2ParamTCneqDimensional.hh"
#include "Framework/ChemicalConsts.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "VariableTransformerSF/VibrEnergy.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Prim2ParamTCneqDimensional, VariableTransformer>
 prim2ParamTcneqDimensionalProv("Prim2ParamTCneqDimensional");

//--------------------------------------------------------------------------//

Prim2ParamTCneqDimensional::Prim2ParamTCneqDimensional(const std::string& objectName) :
 Prim2Param(objectName)
{
}

//--------------------------------------------------------------------------//

Prim2ParamTCneqDimensional::~Prim2ParamTCneqDimensional()
{
}

//--------------------------------------------------------------------------//

void Prim2ParamTCneqDimensional::setup()
{
  LogToScreen(VERBOSE,"Prim2ParamTCneqDimensional::setup() => start\n");

  LogToScreen(VERBOSE,"Prim2ParamTCneqDimensional::setup() => end\n");
}

//--------------------------------------------------------------------------//

void Prim2ParamTCneqDimensional::unsetup()
{
  LogToScreen(VERBOSE,"Prim2ParamTCneqDimensional::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Prim2ParamTCneqDimensional::transform()
{
  LogToScreen(INFO,"Prim2ParamTCneqDimensional::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  double sqrtr;

  // create VibrEnergy object
  VibrEnergy computeVbEnergy;

  rhos.resize( (*nsp) );
  alpha.resize( (*nsp) );
  u.resize(PhysicsInfo::getnbDim());
  T.resize(2);

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {

   for(unsigned ISP=0; ISP<(*nsp); ISP++)
    { rhos.at(ISP) = (*zroe)(ISP,IPOIN); }

   for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++)   
    { u.at(I) = (*zroe)((*nsp)+I,IPOIN); }

   for(unsigned I=0; I<2; I++) 
    { T.at(I) = (*zroe)((*nsp)+PhysicsInfo::getnbDim()+I,IPOIN); }

   // rho
   rho = 0;
   for(unsigned ISP=0; ISP<(*nsp); ISP++) { rho = rho + rhos.at(ISP); }
   sqrtr = sqrt(rho);

   // alpha 
   for(unsigned ISP=0; ISP<(*nsp); ISP++) { alpha.at(ISP) = rhos.at(ISP)/rho; }
   // kinetic energy
   kinetic = pow(u.at(0),2) + pow(u.at(1),2);
   kinetic = kinetic * 0.5;

   // formation enthalpy
   double hftot = 0.0;
   for(unsigned ISP=0; ISP<(*nsp); ISP++) {
    // pow((*uref),2) is used because of the hf dimensionalization
    // made by the ReferenceInfo object
    // in ReferenceInfo hf->at(ISP) = hf->at(ISP)/((*uref) * (*uref));
    double dumhf = hf->at(ISP) * pow(ReferenceInfo::geturef(),2);
    hftot = hftot + alpha.at(ISP) * dumhf; 
   }
   // call for vibrational energy
   computeVbEnergy.callVibrEnergy(T.at(1),alpha);

   ev = computeVbEnergy.getEv();  

   // roto-traslation specif heat
   double Cp = 0;
   for(unsigned ISP=0; ISP<(*nsp); ISP++) {
    Cp = Cp + alpha.at(ISP)* ChemicalConsts::Rgp() / mm->at(ISP) /
         (gams->at(ISP)-1) * gams->at(ISP);
   }

   // total energy
   h = Cp * T.at(0) + ev + hftot + kinetic;

   sqrtr = sqrtr/sqrt(ReferenceInfo::getrhoref());

   (*zroe)((*IEV),IPOIN) = sqrtr * ev /
                (ReferenceInfo::geturef()*ReferenceInfo::geturef());
   (*zroe)((*IX),IPOIN) = sqrtr * u.at(0) / ReferenceInfo::geturef();
   (*zroe)((*IY),IPOIN) = sqrtr * u.at(1) / ReferenceInfo::geturef();
   (*zroe)((*IE),IPOIN) = sqrtr * h / 
               (ReferenceInfo::geturef()*ReferenceInfo::geturef());

   for(unsigned ISP=0; ISP<(*nsp); ISP++) {
    (*zroe)(ISP,IPOIN) = sqrtr * alpha.at(ISP);
   }

   (*XY)(0,IPOIN) = (*XY)(0,IPOIN) / ReferenceInfo::getLref();
   (*XY)(1,IPOIN) = (*XY)(1,IPOIN) / ReferenceInfo::getLref();
  }

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void Prim2ParamTCneqDimensional::transform(vector <double>& m_prim,
                                           vector <double>& m_XY,
                                           vector <double>& m_zroe)
{
  double sqrtr = 0;

  setPhysicsData("FirstCaptured");

  // create VibrEnergy object
  VibrEnergy computeVbEnergy;

  rhos.resize((*nsp),0);
  alpha.resize((*nsp),0);
  u.resize(PhysicsInfo::getnbDim(),0);
  T.resize(2,0);

  for(unsigned ISP=0; ISP<(*nsp); ISP++)
   { rhos.at(ISP) = m_prim.at(ISP); }

  for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++)
   { u.at(I) = m_prim.at((*nsp)+I); }

  for(unsigned I=0; I<2; I++)
   { T.at(I) = m_prim.at((*nsp)+PhysicsInfo::getnbDim()+I); }

  // rho
  rho = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) { rho = rho + rhos.at(ISP); }
  sqrtr = sqrt(rho);

  // alpha 
  for(unsigned ISP=0; ISP<(*nsp); ISP++) { alpha.at(ISP) = rhos.at(ISP)/rho; }

  // kinetic energy
  kinetic = pow(u.at(0),2) + pow(u.at(1),2);
  kinetic = kinetic * 0.5;

  // formation enthalpy 
  double hftot = 0.0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   // pow((*uref),2) is used because of the hf dimensionalization
   // made by the ReferenceInfo object
   // in ReferenceInfo hf->at(ISP) = hf->at(ISP)/((*uref) * (*uref));
   double dumhf = hf->at(ISP) * pow(ReferenceInfo::geturef(),2);
   hftot = hftot + alpha.at(ISP) * dumhf;
  }
  // call for vibrational energy
  computeVbEnergy.callVibrEnergy(T.at(1),alpha);

  ev = computeVbEnergy.getEv();

  // roto-traslation specific heat
  double Cp = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   Cp = Cp + alpha.at(ISP)* ChemicalConsts::Rgp() / mm->at(ISP) /
        (gams->at(ISP)-1) * gams->at(ISP);
  }

  // total energy
  h = Cp * T.at(0) + ev + hftot + kinetic;

  sqrtr = sqrtr/sqrt(ReferenceInfo::getrhoref());

  m_zroe.at((*IEV)) = sqrtr * ev /
               (ReferenceInfo::geturef()*ReferenceInfo::geturef());
  m_zroe.at((*IX)) = sqrtr * u.at(0) / ReferenceInfo::geturef();
  m_zroe.at((*IY)) = sqrtr * u.at(1) / ReferenceInfo::geturef();
  m_zroe.at((*IE)) = sqrtr * h /
              (ReferenceInfo::geturef()*ReferenceInfo::geturef());

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   m_zroe.at(ISP) = sqrtr * alpha.at(ISP);
  }

  m_XY.at(0) = m_XY.at(0) / ReferenceInfo::getLref();
  m_XY.at(1) = m_XY.at(1) / ReferenceInfo::getLref();

}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
