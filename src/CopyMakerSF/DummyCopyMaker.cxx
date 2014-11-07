// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/DummyCopyMaker.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyCopyMaker, CopyMaker>
dummyCopyMakerProv("DummyCopyMaker");

//--------------------------------------------------------------------------//

DummyCopyMaker::DummyCopyMaker(const std::string& objectName) :
  CopyMaker(objectName)
{
}

//--------------------------------------------------------------------------//

DummyCopyMaker::~DummyCopyMaker()
{
}

//--------------------------------------------------------------------------//

void DummyCopyMaker::setup()
{
}

//--------------------------------------------------------------------------//

void DummyCopyMaker::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyCopyMaker::copy()
{
  std::cout << "DummyCopyMaker::copy()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
