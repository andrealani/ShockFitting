// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef OptionValidation_hh
#define OptionValidation_hh

#include "Counter.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

template <typename T> 
class OptionValidation : public Counter {
public:
  
  //constructor
  OptionValidation() {}
  
  //destructor
  virtual ~OptionValidation() {}
  
  // validate the option
  virtual bool isValid(const T& var) = 0;
    
};

//----------------------------------------------------------------------------//

} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
