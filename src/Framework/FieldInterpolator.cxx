// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FieldInterpolator.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//
  
FieldInterpolator::FieldInterpolator(const std::string& objectName) :
  BaseShockFitting(objectName)
{
  m_inXyzIDs = vector<unsigned>();
  addOption("InputXyzIDs", &m_inXyzIDs, "IDs for x,y,z in input field"); 
  
  m_outXyzIDs = vector<unsigned>();
  addOption("OutputXyzIDs", &m_outXyzIDs, "IDs for x,y,z in output (interpolated) field");
}
 
//--------------------------------------------------------------------------//

FieldInterpolator::~FieldInterpolator()
{
}

//--------------------------------------------------------------------------//

void FieldInterpolator::configure(OptionMap& cmap, const std::string& prefix)
{ 
  BaseShockFitting::configure(cmap, prefix);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
