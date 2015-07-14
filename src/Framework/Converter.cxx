// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Converter.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Converter::Converter(const std::string& objectName) :
  BaseShockFitting(objectName)
{
  m_inFmt = "";
  addOption("From",&m_inFmt,
            "Input variables format");
  m_outFmt = "";
  addOption("To",&m_outFmt,
            "Output variables format");
  m_modelTransf = "";
  addOption("GasModel",&m_modelTransf,
            "Gas model used to compute variables transformation");
  m_additionalInfo = "";
  addOption("AdditionalInfo",&m_additionalInfo,
            "CFmesh output and input adimensional or dimensional");
}

//--------------------------------------------------------------------------//

Converter::~Converter()
{
}

//--------------------------------------------------------------------------//

void Converter::configure(OptionMap& cmap, const std::string& prefix)
{
  LogToScreen(VERBOSE, "Converter::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);

  LogToScreen(VERBOSE, "Converter::configure() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
