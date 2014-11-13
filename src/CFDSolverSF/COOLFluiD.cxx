// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFDSolverSF/COOLFluiD.hh"
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
ObjectProvider<COOLFluiD, CFDSolver>
 coolFluiDProv("COOLFluiD");

//--------------------------------------------------------------------------//

COOLFluiD::COOLFluiD(const std::string& objectName) :
 CFDSolver(objectName)
{
}

//--------------------------------------------------------------------------//

COOLFluiD::~COOLFluiD()
{
}

//--------------------------------------------------------------------------//

void COOLFluiD::setup()
{
  LogToScreen(VERBOSE,"COOLFluiD::setup() => start\n");

  LogToScreen(VERBOSE,"COOLFluiD::setup() => end\n");
}

//--------------------------------------------------------------------------//

void COOLFluiD::unsetup()
{
  LogToScreen(VERBOSE,"COOLFluiD::unsetup()\n");
}

//--------------------------------------------------------------------------//

void COOLFluiD::call()
{
  LogToScreen(INFO,"COOLFluiD::call()\n");

  command = "module load cf2-2013.9/solver-openmpi > log/coolfluid.log";
  system(command.c_str());

  if(MeshData::getInstance().getnbProcessors()==1) {
   LogToScreen(DEBUG_MIN,"COOLFluiD::running sequential\n");
   command = "coolfluid-solver --scase ./cf00.CFcase > log/coolfluid.log";
  }

  else if (MeshData::getInstance().getnbProcessors()>1) {
   LogToScreen(DEBUG_MIN,"COOLFluiD::running parallel\n");
   command = "mpi run -np " + MeshData::getInstance().getnbProcessors();
   command = command + " coolfluid-solver --scase ./cf00.CFcase > log/coolfluid.log";
  }

  system(command.c_str());

  if(system(command.c_str())!=0) {
   cout << "COOLFluiD::error => COOLFluiD has return an error code\n"; 
   exit(1); }

}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
