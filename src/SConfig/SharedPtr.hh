// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef SharedPtr_hh
#define SharedPtr_hh

#include <vector>
#include <cstddef>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

template <typename T>
class SharedPtr {
public:
  // constructor 
  SharedPtr(T* ptr = NULL) : m_ptr(ptr) 
  {
    increment();
  }
  
  // copy constructor 
  SharedPtr(const SharedPtr& other) : m_ptr(other.m_ptr) 
  {
    increment();
  }
  
  // reset the pointer
  void reset(T* ptr) 
  {
    decrement();    
    m_ptr = ptr;
    increment();
  }
  
  // destructor 
  ~SharedPtr() {decrement();}
  
  // overloading of operator->
  T* operator->() const {return m_ptr;}
  
  // overloading of operator*
  T& operator*() const {return *m_ptr;}
  
  // overloading operator=
  const SharedPtr& operator=(const SharedPtr<T>& other) 
  {
    reset(other.m_ptr);
    return *this;
  }
  
  // get the bald pointer
  T* get() const {return m_ptr;}
  
private:
  
  // increment the pointer
  void increment() 
  {
    if (m_ptr != NULL) ++(*m_ptr);
  }
  
  // decrement the pointer
  void decrement() 
  {
    if (m_ptr != NULL) {
      --(*m_ptr);	
      if (m_ptr->isZero()) delete m_ptr;
    }
  }
  
private:

  // data ptr
  T* m_ptr;
};

//----------------------------------------------------------------------------//

} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
