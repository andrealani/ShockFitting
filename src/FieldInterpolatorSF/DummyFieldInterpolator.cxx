// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "FieldInterpolatorSF/DummyFieldInterpolator.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyFieldInterpolator, FieldInterpolator> 
dummyFieldInterpolatorProv("DummyFieldInterpolator"); 
 
//--------------------------------------------------------------------------//
  
DummyFieldInterpolator::DummyFieldInterpolator(const std::string& objectName) :
  FieldInterpolator(objectName)
{
}
 
//--------------------------------------------------------------------------//

DummyFieldInterpolator::~DummyFieldInterpolator()
{
}

//--------------------------------------------------------------------------//

void DummyFieldInterpolator::setup()
{
}

//--------------------------------------------------------------------------//

void DummyFieldInterpolator::unsetup() 
{ 
} 
 
//--------------------------------------------------------------------------// 

void DummyFieldInterpolator::interpolate(Field* inField, Field* outField)
{
  std::cout << "DummyFieldInterpolator::interpolate()\n";
}
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting
