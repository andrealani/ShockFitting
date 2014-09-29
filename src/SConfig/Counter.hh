// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Counter_hh
#define Counter_hh

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

class Counter {
public:
  
  // overloading the operator++
  void operator++() {m_counter++;}
  
  // overloading the operator--
  void operator--() {m_counter--;}
  
  // tell if the counter is == 0
  bool isZero() const {return m_counter == 0;}
  
protected:
  
  // constructor
  Counter() : m_counter(0) {}
  
  // destructor
  ~Counter() {}
  
private:
  
  // counter
  int m_counter;
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
