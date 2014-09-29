// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Factory_hh
#define Factory_hh

#include "ConfigException.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

template <class PTYPE> class Provider;

// This class stores an archive of available Provider's for 
// a generic polymorphic type PTYPE. Once registered, they can 
// be accessed by name 
// @author Andrea Lani

template <class PTYPE>
class Factory {
public:
  
  // resister the Provider
  void add(Provider<PTYPE>* ptr) throw()
  {
    if (m_providers.count(ptr->getName()) == 0) {
      m_providers[ptr->getName()] = ptr;
    }
    else {
      throw ConfigException
	(std::string("Provider < " + ptr->getName() + " > is already registered"));
    }
  }
  
  // get the provider corresponding to the given key
  typename PTYPE::PROVIDER* getProvider(const std::string& key) throw()
  {
    if (m_providers.count(key) == 0) throw ConfigException
      (std::string("Provider < " + key + " > is not registered"));
    return dynamic_cast<typename PTYPE::PROVIDER*>(m_providers.find(key)->second);
  }
  
  // instance of the singleton Factory
  static Factory<PTYPE>& getInstance() {static Factory<PTYPE> f; return f;}
  
private:
  
  // constructor
  Factory() {}
  
  // destructor
  ~Factory() 
  { 
  }
  
private: // data
  
  // archive of providers with type <name, Provider*>
  std::map<std::string, Provider<PTYPE>*> m_providers;
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
