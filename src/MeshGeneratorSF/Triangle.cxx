// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>
#include "MeshGeneratorSF/Triangle.hh"
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
ObjectProvider<Triangle, MeshGenerator> triangleProv("Triangle");

//--------------------------------------------------------------------------//

Triangle::Triangle(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

Triangle::~Triangle()
{
}

//--------------------------------------------------------------------------//

void Triangle::setup()
{
  LogToScreen(VERBOSE, "Triangle::setup() => start\n");

  LogToScreen(VERBOSE, "Triangle::setup() => end\n");
}

//--------------------------------------------------------------------------//

void Triangle::unsetup()
{
  LogToScreen(VERBOSE, "Triangle::setup()\n");
}

//--------------------------------------------------------------------------//

void Triangle::generate()
{
  LogToScreen(INFO,"TriangleMeshGenerator::generate()\n");

  fname = MeshData::getInstance().getData<vector<string> >("FNAME");

  command = "../../../src/MeshGeneratorSF/TriangleMeshGeneratorExe/triangle -nep "
            + fname->at(0);
  system(command.c_str());
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
