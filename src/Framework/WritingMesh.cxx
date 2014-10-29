// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/WritingMesh.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

WritingMesh::WritingMesh(const std::string& objectName) :
  BaseShockFitting(objectName)
{
}

//--------------------------------------------------------------------------//

WritingMesh::~WritingMesh()
{

}

//--------------------------------------------------------------------------//

void WritingMesh::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "WritingMesh::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "WritingMesh::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
