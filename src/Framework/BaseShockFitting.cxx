// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "BaseShockFitting.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//
  
BaseShockFitting::BaseShockFitting(const std::string& objectName) : 
  Counter(),
  ConfigObject(objectName)
{
/*  m_flagInputVars = vector<bool>();
  addOption("FlagInputVars",&m_flagInputVars,
            "Flag input variables: 1 (active), 0 (inactive)");
  
  m_flagOutputVars = vector<bool>();
  addOption("FlagOutputVars",&m_flagOutputVars,
            "Flag output variables: 1 (active), 0 (inactive)");
*/
}
 
//--------------------------------------------------------------------------//

BaseShockFitting::~BaseShockFitting()
{
}

//--------------------------------------------------------------------------//

void BaseShockFitting::configure(OptionMap& cmap, const std::string& prefix)
{ 
  ConfigObject::configure(cmap, prefix);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
