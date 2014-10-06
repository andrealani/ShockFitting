// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Remeshing.hh"
#include "Framework/Log.hh"
#include "Framework/IOFunctions.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Remeshing::Remeshing(const std::string& objectName) :
  BaseShockFitting(objectName)
{
}

//--------------------------------------------------------------------------//

Remeshing::~Remeshing()
{
}

//--------------------------------------------------------------------------//

void Remeshing::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "Remeshing::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "Remeshing::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
