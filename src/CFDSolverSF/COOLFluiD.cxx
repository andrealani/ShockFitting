// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFDSolverSF/COOLFluiD.hh"
#include "CFDSolverSF/OverwriteInputFile.hh"
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
  // check if some values of the CoolFluid input file are asked
  // to be changed
  if(m_alterCFDinputfile) {

   stringstream iterValue;
   stringstream oldIterValue;
   stringstream maxNbStepsValue;
   stringstream oldMaxNbStepsValue;

   iterValue.str(string());
   oldIterValue.str(string());
   maxNbStepsValue.str(string());
   oldMaxNbStepsValue.str(string());

   // make a back up of the coolfluid input file
   command.str(string());
   command << "mv cf00.CFcase cf00.CFcase.BAK";
   system(command.str().c_str());
   
   // create the object overwriting the coolfluid file
   OverwriteInputFile ModifyInputCase(string("cf00.CFcase.BAK").c_str(),
                                      string("cf00.CFcase").c_str() );

   // command the object modifying the current coolfluid iter and the 
   // maximum number of coolfluid steps 
   iterValue << "Simulator.SubSystem.InitialIter = " 
             << MeshData::getInstance().getIstep() << "\n";
   oldIterValue << "Simulator.SubSystem.InitialIter = "
             << MeshData::getInstance().getIstep()-1; 
   maxNbStepsValue << "Simulator.SubSystem.MaxNumberSteps.nbSteps = "
             << MeshData::getInstance().getIstep()+1 << "\n";
   oldMaxNbStepsValue << "Simulator.SubSystem.MaxNumberSteps.nbSteps = "
             << MeshData::getInstance().getIstep();
   ModifyInputCase.overwriteValue(oldIterValue.str(),
                                  iterValue.str(),
                                  oldMaxNbStepsValue.str(),
                                  maxNbStepsValue.str());
  }

  LogToScreen(INFO,"COOLFluiD::call()\n");

  if(MeshData::getInstance().getnbProcessors()==1) {
   LogToScreen(DEBUG_MIN,"COOLFluiD::running sequential\n");
   command.str(string());
   command << "coolfluid-solver --scase ./cf00.CFcase >& log/coolfluid.log";
  }

  else if (MeshData::getInstance().getnbProcessors()>1) {
   LogToScreen(DEBUG_MIN,"COOLFluiD::running parallel\n");
   command.str(string());
   command << "  mpirun -np " << MeshData::getInstance().getnbProcessors() ;
   command << " coolfluid-solver --scase ./cf00.CFcase >& log/coolfluid.log";
  }

  FILE* file = popen(command.str().c_str(),"r");

  if(!file) {
   cout << "CFDSolver::error => Failed to call COOLFluiD\n";
   exit(1); }

  pclose(file);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
