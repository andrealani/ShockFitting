// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFDSolverSF/DummyCFDSolver.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyCFDSolver, CFDSolver>
dummyCFDSolverProv("DummyCFDSolver");

//--------------------------------------------------------------------------//

DummyCFDSolver::DummyCFDSolver(const std::string& objectName) :
  CFDSolver(objectName)
{
}

//--------------------------------------------------------------------------//

DummyCFDSolver::~DummyCFDSolver()
{
}

//--------------------------------------------------------------------------//

void DummyCFDSolver::setup()
{
}

//--------------------------------------------------------------------------//

void DummyCFDSolver::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyCFDSolver::call()
{
  std::cout <<"DummyCFDSolver::call()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
