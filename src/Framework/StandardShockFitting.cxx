// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/StandardShockFitting.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/ConfigFileReader.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<StandardShockFitting, ShockFittingObj>
standardShockFittingProv("StandardShockFitting");

//--------------------------------------------------------------------------//

StandardShockFitting::StandardShockFitting(const std::string& objectName) :
  ShockFittingObj(objectName),
  m_readInputFile1(),
  m_readInputFile2()
//  m_BndryNodePtr(),
//  m_RedistrShockPoints()
 // m_FindPhantPoints()
{
}

//--------------------------------------------------------------------------//

StandardShockFitting::~StandardShockFitting()
{
}

//--------------------------------------------------------------------------//

void StandardShockFitting::configure(SConfig::OptionMap& cmap,
                                     const std::string& prefix)
{
  LogToScreen(VERBOSE, "StandardShockFitting::configure() => start\n");

  ShockFittingObj::configure(cmap, prefix);

  LogToScreen(VERBOSE, "StandardShockFitting::configure() => start\n");
}

//--------------------------------------------------------------------------//

void StandardShockFitting::setup()
{
  LogToScreen(VERBOSE, "StandardShockFitting::setup() => start\n");

  ShockFittingObj::setup();

  validate(m_mGenerator.size() == 2,
       "StandardShockFitting::setup() => MeshGeneratorList should have size==2");

//  validate(m_fRemeshing.size() == 2,
     //      "StandardShockFitting::setup() => RemeshingList should have size==3");

  m_readInputFile1 = m_mGenerator[0].ptr();
  m_readInputFile2 = m_mGenerator[1].ptr();
//  m_BndryNodePtr = m_fRemeshing[0].ptr();
//  m_RedistrShockPoints = m_fRemeshing[1].ptr();
 // m_FindPhantPoints = m_fRemeshing[2].ptr();

  LogToScreen(VERBOSE, "StandardShockFitting::setup() => end\n");  
}

//--------------------------------------------------------------------------//

void StandardShockFitting::unsetup()
{
  LogToScreen(VERBOSE, "StandardShockFitting::unsetup() => start\n");

  ShockFittingObj::unsetup();
  
  LogToScreen(VERBOSE, "StandardShockFitting::unsetup() => end\n");
}

//--------------------------------------------------------------------------//

void StandardShockFitting::process()
{
  LogToScreen(VERBOSE, "StandardShockFitting::process() => start\n");

  PhysicsData::getInstance().getPhysicsInfo()->read();
  PhysicsData::getInstance().getChemicalInfo()->read(); 
  PhysicsData::getInstance().getReferenceInfo()->read();

  m_readInputFile1->generate();
  m_readInputFile2->generate();
//  m_BndryNodePtr->remesh();
//  m_RedistrShockPoints->remesh();
//  m_FindPhantPoints->remesh();

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
