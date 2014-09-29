// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/DummyVariableTransformer.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyVariableTransformer, VariableTransformer> 
dummyVariableTransformerProv("DummyVariableTransformer"); 
 
//--------------------------------------------------------------------------//
  
DummyVariableTransformer::DummyVariableTransformer(const std::string& objectName) :
  VariableTransformer(objectName)
{
}
 
//--------------------------------------------------------------------------//

DummyVariableTransformer::~DummyVariableTransformer()
{
}

//--------------------------------------------------------------------------//

void DummyVariableTransformer::setup()
{
}

//--------------------------------------------------------------------------//

void DummyVariableTransformer::unsetup() 
{ 
} 
 
//--------------------------------------------------------------------------// 

void DummyVariableTransformer::transform (double* inArray, double* outArray, int ns)
{
  for (int i = 0; i < ns; ++i) {
    inArray[i] = outArray[i];
  }
  
  std::cout << "DummyVariableTransformer::transform()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
