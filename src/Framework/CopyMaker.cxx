// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CopyMaker.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

CopyMaker::CopyMaker(const std::string& objectName) :
  Counter(),
  ConfigObject(objectName)
{
}

//--------------------------------------------------------------------------//

CopyMaker::~CopyMaker()
{
}

//--------------------------------------------------------------------------//

void CopyMaker::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "CopyMaker::configure() => start\n");

  ConfigObject::configure(cmap, prefix);

  LogToScreen(VERBOSE, "CopyMaker::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

