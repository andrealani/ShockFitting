// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ReferenceInfo.hh"
#include "Framework/ChemicalConsts.hh"
#include "Framework/Log.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

ReferenceInfo::ReferenceInfo(const std::string& objectName) :
  Counter(),
  ConfigObject(objectName)
{
  m_gam = 1;
  addOption("gamma",&m_gam,
            "Isoentropic coefficient of the gas");
  m_Rgas = 1;
  addOption("Rgas",&m_Rgas,
            "Gas constant");
  m_Tref = 1;
  addOption("TempRef",&m_Tref,
            "Reference temperature [K]");
  m_pref = 1;
  addOption("PressRef",&m_pref,
            "Reference pressure [p]");
  m_uref = 1;
  addOption("VelocityRef",&m_uref,
            "Reference speed [m/s]");
  m_rhor = vector<double>();
  addOption("SpeciesDensities",&m_rhor,
            "Species densities");
  m_Lref = 1;
  addOption("Lref",&m_Lref,
            "Reference length");
}

//--------------------------------------------------------------------------//

vector<double> ReferenceInfo::m_rhor;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_Tref=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_pref=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_uref=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_Lref=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_rhoref=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_gam=1;

//--------------------------------------------------------------------------//

double ReferenceInfo::m_Rgas=1;

//--------------------------------------------------------------------------//

ReferenceInfo::~ReferenceInfo()
{
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setup()
{
  LogToScreen(VERBOSE, "ReferenceInfo::setup() => start\n");

  LogToScreen(VERBOSE, "ReferenceInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReferenceInfo::unsetup()
{
  LogToScreen(VERBOSE, "ReferenceInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReferenceInfo::read()
{
  LogToScreen(INFO, "ReferenceInfo::read()\n");

  setPhysicsData();

  logfile.Open(getClassName());

  setReferenceParam();

  if (ChemicalInfo::getModel()=="TCneq") {TCneqAdim();}

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReferenceInfo::TCneqAdim()
{
  Rs->resize(*nsp);
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   hf->at(ISP) = hf->at(ISP)/((m_uref) * (m_uref));
   Rs->at(ISP) = ChemicalConsts::Rgp()/mm->at(ISP);
   Rs->at(ISP) = Rs->at(ISP) * (m_Tref)/((m_uref) * (m_uref));   
  }
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setReferenceParam()
{
  logfile("TREF = ",m_Tref,"\n");
  logfile("PREF = ",m_pref,"\n");
  logfile("UREF = " ,m_uref,"\n");
  logfile("NSP = ",*nsp,"\n");

  if (ChemicalInfo::getModel()=="Cneq" ||
      ChemicalInfo::getModel()== "TCneq") {

   setm_rhoref();

   setAlpha_Rgas_Cv();

   logfile("R = ",R,"\n");
   logfile("Cv = ",Cv,"\n");

   setGamRef();

   logfile("GREF = ",*gref,"\n");

   setGm1Ref();
  }
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setm_rhoref()
{
  m_rhoref = 0;
  assert(m_rhor.size()==(*nsp));
  for (unsigned ISP = 0; ISP<(*nsp); ISP++) {m_rhoref = m_rhoref+m_rhor.at(ISP);}
}

//--------------------------------------------------------------------------//

void ReferenceInfo::setAlpha_Rgas_Cv()
{ 
  double alphaPerc;
  R=0; Cv=0; alpha.resize(*nsp);
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   alpha.at(ISP) = m_rhor.at(ISP)/(m_rhoref);
   R = R + alpha.at(ISP) * ChemicalConsts::Rgp()/mm->at(ISP);
   Cv = Cv + alpha.at(ISP) * ChemicalConsts::Rgp()/mm->at(ISP)/(gams->at(ISP)-1);
   alphaPerc = alpha.at(ISP)*100;
   logfile(namesp->at(ISP),"= ",alphaPerc,"\n");
   logfile("MM = ",mm->at(ISP),"\n");
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
  namesp = PhysicsData::getInstance().getData < vector<string> > ("NAMESP");
  mm = PhysicsData::getInstance().getData < vector<double> > ("MM");
  gams = PhysicsData::getInstance().getData < vector<double> > ("GAMS");
  hf = PhysicsData::getInstance().getData < vector<double> > ("HF");
  Rs = PhysicsData::getInstance().getData < vector<double> > ("RS");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  gm1ref = PhysicsData::getInstance().getData <double> ("GM1REF");
}

//--------------------------------------------------------------------------//

void ReferenceInfo::configure(SConfig::OptionMap& cmap,
                             const std::string& prefix)
{
  LogToScreen(VERBOSE, "ReferenceInfo::configure() => start\n");
  SConfig::ConfigObject::configure(cmap, prefix);
  LogToScreen(VERBOSE, "ReferenceInfo::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
