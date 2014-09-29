// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef SConfigPtr_hh
#define SConfigPtr_hh

#include <string>
#include <memory>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

template <typename OBJ>
class SConfigPtr {
public:
  
  // default constructor
  SConfigPtr() : m_name(""), m_ptr(NULL) {}
  
  // destructor
  ~SConfigPtr() {} 
  
  // get the name
  std::string& name() {return m_name;}
  
  // get the pointer
  OBJ* ptr() const {return m_ptr.get();}
  
  // reset the pointer
  void resetPtr(OBJ* ptr) {m_ptr.reset(ptr);}
  
private:
  
  // self registering/configurable object name
  std::string m_name;
  
  // self registering/configurable object pointer
  std::auto_ptr<OBJ> m_ptr;
  
};

//----------------------------------------------------------------------------//

} // namespace SConfig

//----------------------------------------------------------------------------//

#endif
