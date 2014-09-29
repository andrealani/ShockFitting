// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ConfigObject_hh
#define ConfigObject_hh

#include <vector>

#include "NamedObject.hh"
#include "ConfigMap.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

class ConfigOption;

// This class implements a self-configurable object
// @author Andrea Lani

class ConfigObject : public NamedObject {
public: 
  
  typedef const std::string& ARGS;
  
  // constructor
  ConfigObject(const std::string& name);
  
  // virtual destructor
  virtual ~ConfigObject();
  
  // set the parameter corresponding to the following 
  // @param name       option name in the form "Example.Option. ..."
  // @param var        pointer to the variable to configure
  // @param desc       brief description of the variable functionality
  // @param isDynamic  tells if the option can be configured interactively
  // @param condition  is an option validation object     
  template <typename TYPE>
  void addOption(const std::string& name, TYPE* var, 
		 const std::string& desc, bool isDynamic = false,
		 OptionValidation<TYPE>* condition = NULL)
  {
    const std::string configName1 = getName() + "." + name;
    const std::string configName2 = getParentName() + "." + name;
    ConfigMap::getInstance().add<TYPE>
      (name, configName1, configName2, desc, isDynamic, condition, var, m_options);
  }
  
  // configure all parameters from file
  virtual void configure(OptionMap& cmap, const std::string& prefix);
  
  // configure dependent objects
  virtual void configureDeps(OptionMap& cmap, ConfigObject *const other);
  
protected:
  
  // get the name of the parent
  virtual std::string getParentName() const = 0;
    
private:
  
  // list of options for the current class
  std::vector<ConfigOption*> m_options;
  
  // prefix to prepend to the key
  std::string m_prefix;
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
