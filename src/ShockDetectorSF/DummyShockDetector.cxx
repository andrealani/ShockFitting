// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/DummyShockDetector.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyShockDetector, ShockDetector>
dummyShockDetectorProv("DummyShockDetector");

//--------------------------------------------------------------------------//

DummyShockDetector::DummyShockDetector(const std::string& objectName) :
  ShockDetector(objectName)
{
}

//--------------------------------------------------------------------------//

DummyShockDetector::~DummyShockDetector()
{
}

//--------------------------------------------------------------------------//

void DummyShockDetector::setup()
{
}

//--------------------------------------------------------------------------//

void DummyShockDetector::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyShockDetector::detect()
{
  std::cout << "DummyShockDetector::detect()\n";
}

//--------------------------------------------------------------------------//

void DummyShockDetector::detect(std::vector<double>& firstInfoVector)
{
  std::cout << "DummyShockDetector::detect()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
