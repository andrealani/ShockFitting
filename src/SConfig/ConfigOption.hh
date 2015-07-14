// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ConfigOption_hh
#define ConfigOption_hh 

#include "NamedObject.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements a configurable option
// @author Andrea Lani

class ConfigOption : public NamedObject {
public:
  
  // constructor
  ConfigOption(const std::string& name, bool dyn,
	       const std::string& description) :
    NamedObject(name), 
    m_isDynamic(dyn), m_desc(description)
  {
  }
  
  // destructor
  virtual ~ConfigOption() {}
  
  // configure the value corresponding to this option
  virtual void resetValue(const std::string& value) = 0;
  
  // get the double or int data
  virtual void getData(double* realData, int* intData) = 0;
  
  // description of the option
  std::string description() const {return m_desc;}
  
  // tell if the option is dynamic
  bool isDynamic() const {return m_isDynamic;}
  
  // validate the option
  virtual bool isValid() = 0;
  
  // print value
  virtual void print(std::ofstream& out) const = 0;
  
private:
  
  // flag telling if the option is interactively changeable
  bool m_isDynamic;
  
  // option description
  std::string m_desc;
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
