// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReferenceInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Common/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReferenceInfo, MeshGenerator> referenceInfoProv("ReferenceInfo");

//--------------------------------------------------------------------------//

ReferenceInfo::ReferenceInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{
  m_var = 1;
  addOption("Variables",&m_var,
	    "Variables");
  m_adim = 1;
  addOption("Adimensional",&m_adim,
            "Adimensional");
  m_gam = 1;
  addOption("gamma",&m_gam,
            "Heat Specific ratio");
  m_R = 1;
  addOption("Rgas",&m_R,
            "Gas constant");
  m_Tref = 1;
  addOption("TempRef",&m_Tref,
            "Freestream temperature");
  m_pref = 1;
  addOption("PressRef",&m_pref,
            "Freestream pressure");
  m_uref = 1;
  addOption("VelocityRef",&m_uref,
            "Freestream velocitry");
  m_rhor = vector<double>();
  addOption("SpeciesDensity",&m_rhor,
            "Species density");
  m_Lref = 1;
  addOption("Lref",&m_Lref,
            "L ref");
}

//--------------------------------------------------------------------------//

ReferenceInfo::~ReferenceInfo()
{
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setup()
{
  LogToScreen(VERBOSE, "ReferenceInfo::setup() => start\n");

  setPhysicsData();

  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "ReferenceInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReferenceInfo::unsetup()
{
  LogToScreen(VERBOSE, "ReferenceInfo::unsetup()\n");

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReferenceInfo::generate()
{
  LogToScreen(INFO, "ReferenceInfo::generate()\n");

  setReferenceParam();

  if (model->at(0)=="TCneq") {TCneqAdim();}
}

//--------------------------------------------------------------------------//

void ReferenceInfo::TCneqAdim()
{
  double Rgp = 8.31;
  Rs->resize(*nsp);
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   hf->at(ISP) = hf->at(ISP)/((*uref) * (*uref));
   Rs->at(ISP) = Rgp/mm->at(ISP);
   Rs->at(ISP) = Rs->at(ISP) * (*Tref)/((*uref) * (*uref));   
  }
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setReferenceParam()
{
  *pref = m_pref;
  *Tref = m_Tref;
  *uref = m_uref;
  logfile("TREF = ",*Tref);
  logfile("PREF = ",*pref);
  logfile("UREF = " ,*uref);
  logfile("NSP = ",*nsp);

  if (model->at(0)=="Cneq" || model->at(0) == "TCneq") {

   setRhoRef();

   setAlpha_Rgas_Cv();

   logfile("R = ",R);
   logfile("Cv = ",Cv);

   setGamRef();

   logfile("GREF = ",*gref);

   setGm1Ref();
  }
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setRhoRef()
{
  *rhoref = 0;
  assert(m_rhor.size()==(*nsp));
  for (unsigned ISP = 0; ISP<(*nsp); ISP++) {*rhoref = *rhoref+m_rhor.at(ISP);}
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setAlpha_Rgas_Cv()
{ 
  double Rgp = 8.31;
  double alphaPerc;
  R=0; Cv=0; alpha.resize(*nsp);
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   alpha.at(ISP) = m_rhor.at(ISP)/(*rhoref);
   R = R + alpha.at(ISP) * Rgp/mm->at(ISP);
   Cv = Cv + alpha.at(ISP) * Rgp/mm->at(ISP)/(gams->at(ISP)-1);
   alphaPerc = alpha.at(ISP)*100;
   logfile(namesp->at(ISP),alphaPerc);
   logfile("MM= ",mm->at(ISP));
  }
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setGamRef()
{
  *gref = 1+R/Cv;
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setGm1Ref()
{
  *gm1ref = *gref-1;
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setPhysicsData()
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  model = PhysicsData::getInstance().getData < vector<string> > ("MODEL");
  namesp = PhysicsData::getInstance().getData < vector<string> > ("NAMESP");
  mm = PhysicsData::getInstance().getData < vector<double> > ("MM");
  gams = PhysicsData::getInstance().getData < vector<double> > ("GAMS");
  hf = PhysicsData::getInstance().getData < vector<double> > ("HF");
  Rs = PhysicsData::getInstance().getData < vector<double> > ("RS");
  pref = PhysicsData::getInstance().getData <double> ("PREF");
  Tref = PhysicsData::getInstance().getData <double> ("TREF");
  uref = PhysicsData::getInstance().getData <double> ("UREF");
  rhoref = PhysicsData::getInstance().getData <double> ("RHOREF");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  gm1ref = PhysicsData::getInstance().getData <double> ("GM1REF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
