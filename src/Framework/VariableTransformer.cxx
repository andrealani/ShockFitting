// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/VariableTransformer.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
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
  LogToScreen(VERBOSE, "VariableTransformer::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "VariableTransformer::configure() => end\n");
}

//--------------------------------------------------------------------------//

void VariableTransformer::setup()
{
}

//--------------------------------------------------------------------------//

void VariableTransformer::unsetup() 
{
}

//--------------------------------------------------------------------------//

void VariableTransformer::transform() 
{
}

//--------------------------------------------------------------------------//

void VariableTransformer::transform(std::vector<double>&, std::vector<double>&,
                         std::vector<double>&) 
{
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
