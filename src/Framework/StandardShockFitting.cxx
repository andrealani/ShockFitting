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
  m_readInputFile1(),
  m_readInputFile2(),
  m_readInputFile3(),
  m_readInputValues1(),
  m_readInputValues2(),
  m_readInputValues3(),
  m_remeshField()
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

  validate(m_mGenerator.size() == 6,
           "StandardShockFitting::setup() => MeshGeneratorList should have size==6");

  validate(m_fRemeshing.size() == 1,
           "StandardShockFitting::setup() => RemeshingList should have size==1");

  m_readInputValues1 = m_mGenerator[0].ptr();
  m_readInputValues2 = m_mGenerator[1].ptr();
  m_readInputValues3 = m_mGenerator[2].ptr();
  m_readInputFile1 = m_mGenerator[3].ptr();
  m_readInputFile2 = m_mGenerator[4].ptr();
  m_readInputFile3 = m_mGenerator[5].ptr();
  m_remeshField = m_fRemeshing[0].ptr();

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

  m_readInputValues1->generate();
  m_readInputValues2->generate();
  m_readInputValues3->generate();
  m_readInputFile1->generate();
  m_readInputFile2->generate();
  m_readInputFile3->generate();
  m_remeshField->remesh();

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
