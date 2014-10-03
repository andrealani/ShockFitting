// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PhysicsData.hh"
#include "Common/MeshData.hh"
#include "Framework/StandardShockFitting.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"
#include "MathTools/Array2D.hh"
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
  m_readInputFile(),
  m_readInputValues1(),
  m_readInputValues2()
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

  validate(m_fMeshGenerator.size() == 3,
           "StandardShockFitting::setup() => MeshGeneratorList should have size==2");

  m_readInputFile = m_fMeshGenerator[0].ptr();
  m_readInputValues1 = m_fMeshGenerator[1].ptr();
  m_readInputValues2 = m_fMeshGenerator[2].ptr(); 

  LogToScreen(VERBOSE, "ReadRemeshWrite::setup() => end\n");  
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

  m_readInputFile->generate();
  m_readInputValues1->generate();
  m_readInputValues2->generate();

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
