// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/PhysicsInfo.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

PhysicsInfo::PhysicsInfo(const std::string& objectName) :
  Counter(),
  ConfigObject(objectName)
{
  m_ndim = 1;
  addOption("NDIM",&m_ndim,
            "Space dimension");
  m_gam = 1;
  addOption("GAM",&m_gam,
            "Heat specific ratio");
  m_ndofmax = 1;
  addOption("NDOFMAX",&m_ndofmax,
            "Maximum number of degree of freedom");

  m_nshmax = 1;
  addOption("NSHMAX",&m_nshmax,
            "Maximum number of shocks");

  m_npshmax = 1;
  addOption("NPSHMAX",&m_npshmax,
            "Maximum number of shock points for each shock");

  m_neshmax = 1;
  addOption("NESHMAX",&m_neshmax,
            "Maximum number od shock edges for each shock");

  m_nspmax = 1;
  addOption("NSPMAX",&m_nspmax,
            "Maximum number of special points");

  m_naddholesmax = 1;
  addOption("NADDHOLESMAX",&m_naddholesmax,
            "Maximum number of holes");
}

//--------------------------------------------------------------------------//

PhysicsInfo::~PhysicsInfo()
{
}

//--------------------------------------------------------------------------//

void PhysicsInfo::setup()
{
  LogToScreen(VERBOSE, "PhysicsInfo::setup() => start\n");

  setPhysicsData();

  LogToScreen(VERBOSE, "PhysicsInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void PhysicsInfo::unsetup()
{
  LogToScreen(VERBOSE, "PhysicsInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void PhysicsInfo::read()
{
  LogToScreen(INFO, "PhysicsInfo::read()\n");

  *ndim = m_ndim;
  *gam = m_gam;
  *gm1 = *gam-1;
  *ndofmax = m_ndofmax;
  *npshmax = m_npshmax;
  *nshmax = m_nshmax;
  *neshmax = m_neshmax;
  *nspmax = m_nspmax;
  *naddholesmax = m_naddholesmax;
}

//--------------------------------------------------------------------------//

void PhysicsInfo::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  gam = PhysicsData::getInstance().getData <double> ("GAM");
  gm1 = PhysicsData::getInstance().getData <double> ("GM1");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  neshmax = PhysicsData::getInstance().getData <unsigned> ("NESHMAX");
  nspmax = PhysicsData::getInstance().getData <unsigned> ("NSPMAX");
  naddholesmax = PhysicsData::getInstance().getData <unsigned> ("NADDHOLESMAX");
}

//--------------------------------------------------------------------------//

void PhysicsInfo::configure(SConfig::OptionMap& cmap,
                             const std::string& prefix)
{
  LogToScreen(VERBOSE, "PhysicsInfo::configure() => start\n");
  SConfig::ConfigObject::configure(cmap, prefix);
  LogToScreen(VERBOSE, "PhysicsInfo::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
