// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>
#include "MeshGeneratorSF/TriangleExe.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<TriangleExe, MeshGenerator> triangleExeProv("TriangleExe");

//--------------------------------------------------------------------------//

TriangleExe::TriangleExe(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

TriangleExe::~TriangleExe()
{
}

//--------------------------------------------------------------------------//

void TriangleExe::setup()
{
  LogToScreen(VERBOSE, "TriangleExe::setup() => start\n");

  LogToScreen(VERBOSE, "TriangleExe::setup() => end\n");
}

//--------------------------------------------------------------------------//

void TriangleExe::unsetup()
{
  LogToScreen(VERBOSE, "TriangleExe::unsetup()\n");
}

//--------------------------------------------------------------------------//

void TriangleExe::generate()
{
  LogToScreen(INFO,"TriangleExe::generate()\n");

  fname = MeshData::getInstance().getData <stringstream>("FNAME");
  command = "/data/deamicis/ShockFitting.git/trunk/src/MeshGeneratorSF/TriLibrary/triangle -nep "
            + fname->str() + " > log/TriangleExe.log";

  system(command.c_str());

  if(system(command.c_str())!=0) {
   cout << "TriangleExe::error => Triangle Mesh Generator execution failed\n";
   exit(1); }
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
