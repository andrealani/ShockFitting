// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ObjectProvider_hh
#define ObjectProvider_hh

#include <string>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements a generic provider allowing client code
// to create arbitrary types of polymorphic objects. 
// @author Andrea Lani
template <class CTYPE, class PTYPE>
class ObjectProvider : public PTYPE::PROVIDER {
public:
  
  // constructor
  ObjectProvider(const std::string& name) : PTYPE::PROVIDER(name) {}
  
  // destructor
  ~ObjectProvider() {}
  
  // create an  object of polymorphic type PTYPE (parent) and
  // static type CTYPE (concrete)
  PTYPE* create(typename PTYPE::ARGS args) {return new CTYPE(args);}
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
