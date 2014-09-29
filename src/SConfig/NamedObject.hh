// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef NamedObject_hh
#define NamedObject_hh

#include <string>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class is meant to be inherited from. 
// It provides a name to the deriving object.
// @author Andrea Lani
class NamedObject {
public:
  
  // constructor
  NamedObject(const std::string& name) : m_name(name) {}
  
  // destructor
  ~NamedObject() {}
  
  // get the name
  std::string getName() const {return m_name;}
  
private:
  
  // name of the object
  std::string m_name;
};

//----------------------------------------------------------------------------//

} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
