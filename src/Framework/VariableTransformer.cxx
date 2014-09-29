// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/VariableTransformer.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//
  
VariableTransformer::VariableTransformer(const std::string& objectName) :
  BaseShockFitting(objectName)
{
}
 
//--------------------------------------------------------------------------//

VariableTransformer::~VariableTransformer()
{
}

//--------------------------------------------------------------------------//

void VariableTransformer::configure(OptionMap& cmap, const std::string& prefix)
{ 
  BaseShockFitting::configure(cmap, prefix);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
