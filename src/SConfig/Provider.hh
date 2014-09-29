// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Provider_hh
#define Provider_hh

#include "Factory.hh"
#include "NamedObject.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements an object provider.
// @author Andrea Lani

template <class PTYPE>
class Provider : public NamedObject {
public:
  
  // constructor
  Provider(const std::string& name) : NamedObject(name)
  {
    std::cout << "Registering object named [" << name << "]\n"; 
    Factory<PTYPE>::getInstance().add(this);
  }
  
  // destructor
  virtual ~Provider() {}
  
  // Create the chosen concrete object
  virtual PTYPE* create(typename PTYPE::ARGS args) = 0;
  

};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
