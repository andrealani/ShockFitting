// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConfigObject.hh"
#include "ConfigOption.hh"
#include "ConfigFileReader.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

ConfigObject::ConfigObject(const std::string& name) :
  NamedObject(name), 
  m_options(),
  m_prefix()
{
}

//----------------------------------------------------------------------------//

ConfigObject::~ConfigObject()
{
}

//----------------------------------------------------------------------------//

// configure all parameters from file
void ConfigObject::configure(OptionMap& cmap, const string& prefix)
{  
  const bool first = ConfigFileReader::isFirstConfig();
  
  // back the prefix up
  m_prefix = prefix;
  // cout << "prefix = " << prefix << endl;
  
  for (unsigned i = 0; i < m_options.size(); ++i) {
    if (first || (!first && m_options[i]->isDynamic())) {
      // cout << m_options[i]->getName() << ", dynamic ? " <<  m_options[i]->isDynamic() << endl;
      const string key = prefix + m_options[i]->getName();
      // cout << "key = " << key << endl;
      
      bool found = false;
      pair<OptionMap::MapItr, OptionMap::MapItr> pIt = cmap.find(key, found);
      if (found) { 
	pair<string, bool*> p = pIt.first->second;
	m_options[i]->resetValue(p.first);
	
	// set the boolean flag to true, so that we know that this option 
	// has been actually used
	assert(p.second != NULL);
	(*p.second) = true; 
      }
    }
  }
}
  
//----------------------------------------------------------------------------//

// configure dependent objects
void ConfigObject::configureDeps(OptionMap& cmap, ConfigObject* other) 
{
  const string prefix = m_prefix + getName() + "."; 
  other->configure(cmap, prefix);
}

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//
