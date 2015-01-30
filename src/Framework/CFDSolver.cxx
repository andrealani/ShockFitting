// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CFDSolver.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

CFDSolver::CFDSolver(const std::string& objectName) :
  BaseShockFitting(objectName)
{
  m_alterCFDinputfile = "false";
  addOption("alterCFDinputFile",&m_alterCFDinputfile,
            "Ask to change some values in the CFD input file at run time");
  m_alterWhichValues = vector<string>(1);
  addOption("alterWhichValues",&m_alterWhichValues,
            "Specify the variables to be changed");
}

//--------------------------------------------------------------------------//

CFDSolver::~CFDSolver()
{
}

//--------------------------------------------------------------------------//

void CFDSolver::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "CFDSolver::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "CFDSolver::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

