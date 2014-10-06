// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <vector>
#include "MeshGeneratorSF/PhysicsInfo.hh"
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
ObjectProvider<PhysicsInfo, MeshGenerator> physicsInfoProv("PhysicsInfo");

//--------------------------------------------------------------------------//

PhysicsInfo::PhysicsInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{

  m_ndim = 1;
  addOption("NDIM",&m_ndim,
            "Space dimension");

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

void PhysicsInfo::generate()
{
  LogToScreen(INFO, "PhysicsInfo::generate()\n");

  *ndim = m_ndim;
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
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  neshmax = PhysicsData::getInstance().getData <unsigned> ("NESHMAX");
  nspmax = PhysicsData::getInstance().getData <unsigned> ("NSPMAX");
  naddholesmax = PhysicsData::getInstance().getData <unsigned> ("NADDHOLESMAX");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
