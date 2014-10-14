// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/DummyMeshGenerator.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyMeshGenerator, MeshGenerator>
dummyMeshGeneratorProv("DummyMeshGenerator");

//--------------------------------------------------------------------------//

DummyMeshGenerator::DummyMeshGenerator(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

DummyMeshGenerator::~DummyMeshGenerator()
{
}

//--------------------------------------------------------------------------//

void DummyMeshGenerator::setup()
{
}

//--------------------------------------------------------------------------//

void DummyMeshGenerator::unsetup()
{
}

//--------------------------------------------------------------------------// 

void DummyMeshGenerator::generate()
{
  std::cout << "DummyMeshGenerator::generate()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
