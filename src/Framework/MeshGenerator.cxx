// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MeshGenerator.hh"
#include "Common/PhysicsData.hh"
#include "Framework/Log.hh"
#include "Framework/IOFunctions.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

MeshGenerator::MeshGenerator(const std::string& objectName) :
  BaseShockFitting(objectName)
{
  m_inputFiles = vector<string>();
  addOption("InputFiles",&m_inputFiles,
            "List of the names of input files");

}


//--------------------------------------------------------------------------//

MeshGenerator::~MeshGenerator()
{
}

//--------------------------------------------------------------------------//

void MeshGenerator::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "MeshGenerator::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "MeshGenerator::configure() => end\n");
}

//--------------------------------------------------------------------------//

void MeshGenerator::setMeshField(Field *const field)
{
  LogToScreen(VERBOSE, "MeshGenerator::setMeshField()\n");
}

//--------------------------------------------------------------------------//

void MeshGenerator::getMeshField(Field* field)
{
  LogToScreen(VERBOSE, "MeshGenerator::getMeshField()\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
