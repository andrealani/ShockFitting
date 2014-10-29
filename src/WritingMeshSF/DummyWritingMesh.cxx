// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "WritingMeshSF/DummyWritingMesh.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyWritingMesh, WritingMesh>
dummyWritingMeshProv("DummyWritingMesh");

//--------------------------------------------------------------------------//

DummyWritingMesh::DummyWritingMesh(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

DummyWritingMesh::~DummyWritingMesh()
{
}

//--------------------------------------------------------------------------//

void DummyWritingMesh::setup()
{
}

//--------------------------------------------------------------------------//

void DummyWritingMesh::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyWritingMesh::write()
{
  std::cout << "DummyWritingMesh::write()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

