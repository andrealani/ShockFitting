// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef StringManip_hh
#define StringManip_hh

#include <string>
#include <sstream>
#include <map>

#include "ConfigException.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

template <typename T>
static std::string to_str(T v)
{
  std::ostringstream oss; oss << v; return oss.str();
}

//----------------------------------------------------------------------------//

template <typename T>
static T fromString (const std::string& s)
{
  if (s.size() > 0) {
    T v;
    std::istringstream iss(s.c_str()); iss >> v;
    return v;
  }
  return T();
}

//----------------------------------------------------------------------------//

// trim front/rear of a string
static void trim (std::string& s)
{
  const size_t start = s.find_first_not_of(" \t");
  const size_t end = s.find_last_not_of(" \t");
  if ((start == std::string::npos) || (end == std::string::npos)) {s = "";}
  else {s = s.substr(start, end - start+1);}
}

//----------------------------------------------------------------------------//

// trim front of a string
static void trimFront (std::string& s)
{
  const size_t pos = s.find_first_not_of(" \t");
  if (pos != std::string::npos) s = s.substr(pos);
}

//----------------------------------------------------------------------------//

// trim rear of a string
static void trimRear (std::string& s)
{
  const size_t pos = s.find_last_not_of(" \t");
  if (pos != std::string::npos) s = s.substr(0, pos+1);
}

//----------------------------------------------------------------------------//

/// Substitute the rhs for the lhs wherever it appears in the string
/// @param lhs partial std::string to match
/// @param rhs partial std::string to substitute
/// @param out string to modify
static void subst (const std::string& lhs, const std::string& rhs, std::string& out)
{ 
  std::string::size_type i = 0;
  while (i < out.length()) {
    std::string substr (out, i);
    std::string::size_type idx = out.find(lhs, i);
    if (idx == std::string::npos) {
      break;  // we're at the end 
    }
    else {
      out.replace(idx, lhs.length(), rhs);
      i = idx + rhs.length();
      if (!rhs.length())
	++i;
    }
  }
}

//----------------------------------------------------------------------------//                                                                                                              

// Proxy for a std::map for StringInt
class StringIntMap : public std::map<std::string, int> {
public:  
  
  // static function making this calss a singleton
  static StringIntMap& getInstance() 
  {
    static StringIntMap s;
    return s;    
  }
  
  // find key
  int find(const std::string& s) 
  {
    if (this->count(s) == 0) throw ConfigException
      (std::string("StringIntMap < " + s + " > is not available"));
    return std::map<std::string, int>::find(s)->second; 
  }
    
private:
  
  StringIntMap() : std::map<std::string, int>() {}
  ~StringIntMap() {}
};

//----------------------------------------------------------------------------//

class StringInt {
public:
  
  // accessor/mutator 
  int& oint() {return m_oint;}
  
  // accessor/mutator 
  int oint() const {return m_oint;}
  
  // accessor/mutator
  std::string& name() {return m_name;}
  
  // accessor/mutator
  std::string name() const {return m_name;}
  
  // overloading of assignment operator for input StringInt
  const StringInt& operator=(const StringInt& s) 
  {
    m_oint = s.m_oint;
    m_name = s.m_name;
    return *this;
  }
  
  // overloading of assignment operator for input string
  const StringInt& operator=(const std::string& s) 
  {    
    m_oint = StringIntMap::getInstance().find(s);
    m_name = s;
    return *this;
  }
  
  // overloading of assignment operator for input int
  const StringInt& operator=(int oint) 
  {    
    m_oint = oint;
    m_name = to_str(oint);
    return *this;
  }
  
  // overloading of the output operator
  friend std::ostream& operator<< (std::ostream& out, const StringInt& s) 
  {
    out << "{" << s.m_oint << " - " << s.m_name << "}\n";
    return out;
  }
  
  // conversion operator
  operator int () {return m_oint;}
  
  // overloading of the operators ==, !=, >=, <=, <, >
#define STRINGINT_OP_INT(__op__) \
friend bool operator __op__(const StringInt& s, int oin) \
{return (s.m_oint __op__ oin);}
  
STRINGINT_OP_INT(==)
STRINGINT_OP_INT(!=)
STRINGINT_OP_INT(>=)
STRINGINT_OP_INT(<=)
STRINGINT_OP_INT(>)
STRINGINT_OP_INT(<)

#undef STRINGINT_OP_INT
  
private:

  std::string m_name;
  int m_oint;
  
};

//----------------------------------------------------------------------------//

template <typename T>
class StringT {
public:
  
  // accessor/mutator 
  T& ptr() {return m_ptr;}
  
  // accessor/mutator 
  T ptr() const {return m_ptr;}
  
  // accessor/mutator
  std::string& name() {return m_name;}
  
  // accessor/mutator
  std::string name() const {return m_name;}
  
  // overloading of assignment operator for input StringT
  const StringT& operator=(const StringT& s) 
  {
    m_ptr = s.m_ptr;
    m_name = s.m_name;
    return *this;
  }
  
  // overloading of assignment operator for input string
  const StringT& operator=(const std::string& s) 
  {
    m_name = s;
    return *this;
  }
  
  // overloading of the output operator
  friend std::ostream& operator<< (std::ostream& out, const StringT& s) 
  {
    out << "{ name = " << s.m_name << "}\n";
    return out;
  }
  
// overloading of the operators ==, !=, >=, <=, <, >
#define STRINGT_OP_INT(__op__) \
friend bool operator __op__(const StringT& s1, const StringT& s2) \
{return (s1.m_name __op__ s2.m_name);}
  
STRINGT_OP_INT(==)
STRINGT_OP_INT(!=)
STRINGT_OP_INT(>=)
STRINGT_OP_INT(<=)
STRINGT_OP_INT(>)
STRINGT_OP_INT(<)

#undef STRINGT_OP_INT
  
private:

  std::string m_name;
  T m_ptr;
  
};

//----------------------------------------------------------------------------//

} // end namespace SConfig 

//----------------------------------------------------------------------------//

#endif
