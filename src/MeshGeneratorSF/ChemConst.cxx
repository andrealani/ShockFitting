// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ChemConst.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Common/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ChemConst, MeshGenerator> chemConstProv("ChemConst");

//--------------------------------------------------------------------------//

ChemConst::ChemConst(const std::string& objectName) :
  MeshGenerator(objectName)
{
  m_R = 1;
  addOption("R",&m_R,
            "Gas constant");
 
  m_Na = 1;
  addOption("Na",&m_Na,
	    "Avogadro constant");

  m_K = 1;
  addOption("K",&m_K,
	    "Boltzmann constant");

  m_KSI = 1;
  addOption("KSI",&m_KSI,
	    "Boltzmann constant (International System");
}

//--------------------------------------------------------------------------//

ChemConst::~ChemConst()
{
}

//--------------------------------------------------------------------------//

void ChemConst::setup()
{
  LogToScreen(VERBOSE, "ChemConst::setup() => start\n");

  setPhysicsData();

  LogToScreen(VERBOSE, "ChemConst::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ChemConst::unsetup()
{
  LogToScreen(VERBOSE, "ChemConst::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ChemConst::generate()
{
  LogToScreen(INFO, "ChemConst::generate()\n");

  *r_R = m_R;
  *r_Na = m_Na;
  *r_K = m_K;
  *r_KSI = m_KSI;
}

//--------------------------------------------------------------------------//

void ChemConst::setPhysicsData()
{
  r_R = PhysicsData::getInstance().getData <double> ("R");
  r_Na = PhysicsData::getInstance().getData <double> ("Na");
  r_K = PhysicsData::getInstance().getData <double> ("K");
  r_KSI = PhysicsData::getInstance().getData <double> ("KSI");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
