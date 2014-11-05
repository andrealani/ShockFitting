// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/DummyStateUpdater.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyStateUpdater, StateUpdater>
dummyStateUpdaterProv("DummyStateUpdater");

//--------------------------------------------------------------------------//

DummyStateUpdater::DummyStateUpdater(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//--------------------------------------------------------------------------//

DummyStateUpdater::~DummyStateUpdater()
{
}

//--------------------------------------------------------------------------//

void DummyStateUpdater::setup()
{
}

//--------------------------------------------------------------------------//

void DummyStateUpdater::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyStateUpdater::update()
{
  std::cout << "DummyStateUpdater::update\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
