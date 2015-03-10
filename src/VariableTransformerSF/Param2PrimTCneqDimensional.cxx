// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Param2PrimTCneqDimensional.hh"
#include "Framework/ChemicalConsts.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/ReferenceInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "VariableTransformerSF/ComputeTv.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Param2PrimTCneqDimensional, VariableTransformer>
 param2PrimTCneqDimensionalProv("Param2PrimTCneqDimensional");

//--------------------------------------------------------------------------//

Param2PrimTCneqDimensional::Param2PrimTCneqDimensional
                            (const std::string& objectName) :
 Param2Prim(objectName)
{
}

//--------------------------------------------------------------------------//

Param2PrimTCneqDimensional::~Param2PrimTCneqDimensional()
{
}

//--------------------------------------------------------------------------//

void Param2PrimTCneqDimensional::setup()
{
  LogToScreen(VERBOSE,"Param2PrimTCneqDimensional::setup() => start\n");

  LogToScreen(VERBOSE,"Param2PrimTCneqDimensional::setup() => end\n");
}

//--------------------------------------------------------------------------//

void Param2PrimTCneqDimensional::unsetup()
{
  LogToScreen(VERBOSE,"Param2PrimTCneqDimensional::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Param2PrimTCneqDimensional::transform()
{
  LogToScreen(INFO,"Param2PrimTCneqDimensional::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    // zrho and rho
    double sqrtr = 0;
    for (unsigned ISP=0; ISP<(*nsp); ISP++) {
     sqrtr = sqrtr + (*zroe)(ISP,IPOIN);  }
    rho = pow(sqrtr,2);

    // rho-i and alpha-i
    alpha.resize(*nsp);
    rhos.resize(*nsp);
    for (unsigned ISP=0; ISP<(*nsp); ISP++) {
     alpha.at(ISP) = (*zroe)(ISP,IPOIN)/sqrtr;
     rhos.at(ISP) = (*zroe)(ISP,IPOIN) * sqrtr * ReferenceInfo::getrhoref();
 }     

    // u, v, h, ev
    u.resize(2);
    h = (*zroe)((*IE),IPOIN)/sqrtr * pow(ReferenceInfo::geturef(),2);
    u.at(0) = (*zroe)((*IX),IPOIN)/sqrtr * ReferenceInfo::geturef();
    u.at(1) = (*zroe)((*IY),IPOIN)/sqrtr * ReferenceInfo::geturef();
    ev = (*zroe)((*IEV),IPOIN)/sqrtr * 
         pow(ReferenceInfo::geturef(),2);

    // kinetic energy
    kinetic = pow(u.at(0),2)+pow(u.at(1),2);
    kinetic = kinetic * 0.5;

    // formation enthalpy
    double hftot = 0;
    for(unsigned ISP=0; ISP<(*nsp); ISP++) {
     // pow((*uref),2) is used because of the hf dimensionalization
     // done by the ReferenceInfo object
     // in ReferenceInfo hf->at(ISP) = hf->at(ISP)/((*uref) * (*uref));
     hftot = hftot + alpha.at(ISP) * hf->at(ISP)*
             pow(ReferenceInfo::geturef(),2);
    }

    // roto-traslation enthalpy
    double htr = h - kinetic - hftot - ev;

    // roto-translation specific heat
    double Cp = 0;
    for(unsigned ISP=0; ISP<(*nsp); ISP++) {
     Cp = Cp + alpha.at(ISP) * gams->at(ISP)/(gams->at(ISP)-1)
          * ChemicalConsts::Rgp()/mm->at(ISP);
    }

    // temperature
    T.resize(2);

    T.at(0) = htr/Cp;

    // vibrational temperature
    ComputeTv cTv;
    if((*nmol)==1) {

     T.at(1) = T.at(0)*pow(10,-6);

     cTv.callComputeTv(ev,alpha,T);
     T = cTv.getT();
    }

    else if (*nmol>1) { T.at(1)=T.at(0);
                        cTv.callComputeTv(ev,alpha,T);
                        T=cTv.getT(); }

    else { cout << "Parm2PrimTCneq::error => NMOL = 0\n";
           exit(1); }

    for(unsigned ISP=0; ISP<(*nsp); ISP++) {
     (*zroe)(ISP,IPOIN)=rhos.at(ISP);
    }

   (*zroe)((*IE),IPOIN) = u.at(0);
   (*zroe)((*IX),IPOIN) = u.at(1);
   (*zroe)((*IY),IPOIN) = T.at(0);
   (*zroe)((*IEV),IPOIN) = T.at(1);

   (*XY)(0,IPOIN) = (*XY)(0,IPOIN) * ReferenceInfo::getLref();
   (*XY)(1,IPOIN) = (*XY)(1,IPOIN) * ReferenceInfo::getLref();
  }

  // de-allocate dynamic arrays
  freeArray();
}  

//--------------------------------------------------------------------------//

void Param2PrimTCneqDimensional::transform(vector <double>* m_zroe,
					   vector <double>* m_XY,
                                           vector <double>* m_prim)
{
  setPhysicsData();

  // zrho and rho
  double sqrtr = 0;
  for (unsigned ISP=0; ISP<(*nsp); ISP++) { 
   sqrtr = sqrtr + m_zroe->at(ISP);  }
   rho = pow(sqrtr,2);
 
  // rho-i and alpha-i
  alpha.resize(*nsp);
  rhos.resize(*nsp);
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   alpha.at(ISP) = m_zroe->at(ISP)/sqrtr;
   rhos.at(ISP) = m_zroe->at(ISP) * sqrtr * ReferenceInfo::getrhoref();
  } 

  // u, v, h, ev
  u.resize(2);
  h = m_zroe->at((*IE))/sqrtr * pow(ReferenceInfo::geturef(),2);
  u.at(0) = m_zroe->at((*IX))/sqrtr * ReferenceInfo::geturef();
  u.at(1) = m_zroe->at((*IY))/sqrtr * ReferenceInfo::geturef();
  ev = m_zroe->at((*IEV))/sqrtr * pow(ReferenceInfo::geturef(),2);

  // kinetic energy
  kinetic = pow(u.at(0),2)+pow(u.at(1),2);
  kinetic = kinetic * 0.5;

  // formation enthalpy
  double hftot = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   // pow((*uref),2) is used because of the hf dimensionalization
   // done by the ReferenceInfo object
   // in ReferenceInfo hf->at(ISP) = hf->at(ISP)/((*uref) * (*uref));
   hftot = hftot + alpha.at(ISP) * hf->at(ISP)*pow(ReferenceInfo::geturef(),2);
  }

  // roto-traslation enthalpy
  double htr = h - kinetic - hftot - ev;

  // roto-translation specific heat
  double Cp = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   Cp = Cp + alpha.at(ISP) * gams->at(ISP)/(gams->at(ISP)-1)
        * ChemicalConsts::Rgp()/mm->at(ISP);
  }

  // temperature
  T.resize(2);

  T.at(0) = htr/Cp;

  // vibrational temperature
  ComputeTv cTv;
  if((*nmol)==1) {

   T.at(1) = T.at(0)*pow(10,-6);

   cTv.callComputeTv(ev,alpha,T);
   T = cTv.getT();
  }

  else if (*nmol>1) { T.at(1)=T.at(0);
                      cTv.callComputeTv(ev,alpha,T);
                      T=cTv.getT(); }

  else { cout << "ConverterSF::Parm2Prim4TCneq::error => NMOL = 0\n";
         exit(1); }

  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   m_prim->at(ISP)=rhos.at(ISP);
  }

  m_prim->at((*IE)) = u.at(0);
  m_prim->at((*IX)) = u.at(1);
  m_prim->at((*IY)) = T.at(0);
  m_prim->at((*IEV)) = T.at(1);

  m_XY->at(0) = m_XY->at(0) * ReferenceInfo::getLref();
  m_XY->at(1) = m_XY->at(1) * ReferenceInfo::getLref();
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

