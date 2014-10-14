// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/DummyRemeshing.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyRemeshing, Remeshing>
dummyRemeshingProv("DummyRemeshing");

//--------------------------------------------------------------------------//

DummyRemeshing::DummyRemeshing(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

DummyRemeshing::~DummyRemeshing()
{
}

//--------------------------------------------------------------------------//

void DummyRemeshing::setup()
{
}

//--------------------------------------------------------------------------//

void DummyRemeshing::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyRemeshing::remesh()
{
  std::cout << "DummyRemeshing::remesh()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
