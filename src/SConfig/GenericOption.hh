// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef GenericOption_hh
#define GenericOption_hh 

#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <vector>

#include "StringManip.hh"
#include "ConfigOption.hh"
#include "OptionValidation.hh"
#include "SharedPtr.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements a template computational 
// constructor that assigns data to one of the two arguments 
// double* or int*, according to the type.
// Different specializations are available.
// @author Andrea Lani

template <typename T>
struct GetDataFun 
{
  // default behaviour is doing nothing
  GetDataFun(double* realData, int* intData, T* var) {}
};

template <>
struct GetDataFun<double> {
  GetDataFun(double* realData, int* intData, double* var) 
  {assert(var != NULL); assert(realData != NULL); *realData = *var;}
};

template <>
struct GetDataFun<int> {
  GetDataFun(double* realData, int* intData, int* var) 
  {assert(var != NULL); assert(intData != NULL); *intData = *var;}
};

template <>
struct GetDataFun<StringInt> {
  GetDataFun(double* realData, int* intData, StringInt* var) 
  {assert(var != NULL); assert(intData != NULL); *intData = var->oint();}
};

template <>
struct GetDataFun< std::vector<double> > {
  GetDataFun(double* realData, int* intData, std::vector<double>* var) 
  {
    assert(var != NULL);  assert(realData != NULL);
    std::copy(&(*var)[0], &(*var)[0]+var->size(), &realData[0]);
  }
};

template <>
struct GetDataFun< std::vector<int> > {
  GetDataFun(double* realData, int* intData, std::vector<int>* var) 
  {
    assert(var != NULL); assert(intData != NULL);
    std::copy(&(*var)[0], &(*var)[0]+var->size(), &intData[0]);
  }
};

template <>
struct GetDataFun< std::vector<StringInt> > {
  GetDataFun(double* realData, int* intData, std::vector<StringInt>* var) 
  {
    assert(var != NULL); assert(intData != NULL);
    for (unsigned i = 0; i < var->size(); ++i) {
      intData[i] = (*var)[i].oint();
    }
  }
};

//----------------------------------------------------------------------------//

template <typename T>
struct ConvertFromString {
  // default behaviour is doing nothing
  ConvertFromString(T* var, const std::string& value) 
  {
    *var = fromString<T>(value);
  }
};

template <>
struct ConvertFromString<StringInt> {
  // default behaviour is doing nothing
  ConvertFromString(StringInt* var, const std::string& value) 
  {
    var->oint() = StringIntMap::getInstance().find(value);
    var->name() = value;
  }
};

template <typename T>
struct ConvertFromString< StringT<T> > {
  // default behaviour is doing nothing
  ConvertFromString(StringT<T>* var, const std::string& value) 
  {
    var->ptr()  = T();
    var->name() = value;
  }
};

//----------------------------------------------------------------------------//

// This class implements a generic option to be configured automatically
// @author Andrea Lani

template <typename T>
class GenericOption : public ConfigOption {
public:
  
  // constructor
  GenericOption(const std::string& name, bool dyn, 
		const std::string& description, 
		OptionValidation<T>* condition, T* var) : 
    ConfigOption(name, dyn, description), 
    m_var(var),
    m_condition(condition)
  {
  }
  
  // destructor
  ~GenericOption() {}
  
  // get the double or int data
  void getData(double* realData, int* intData) 
  {
    GetDataFun<T>(realData, intData, m_var);
  }
  
  // configure the value corresponding to this option
  void resetValue(const std::string& value) 
  { 
    // convert string to value
    ConvertFromString<T>(m_var, value);
  }
  
  // validate the option
  bool isValid() 
  {
    if (m_condition.get() != NULL) {
      assert(m_var != NULL);
      return m_condition->isValid(*m_var);
    }
    return true;
  }
  
  // print value
  void print(std::ofstream& out) const 
  {
    out << *m_var << std::endl;
  }
  
private:
  
  // variable to configure
  T* m_var;
  
  // option validity condition
  SharedPtr< OptionValidation<T> > m_condition;
  
};

//----------------------------------------------------------------------------//

// This class implements a partial specialization for std::vector types
// for a generic option to be configured automatically
// @author Andrea Lani

template <typename T>
class GenericOption< std::vector<T> > : public ConfigOption {
public:
  
  // constructor
  GenericOption(const std::string& name, bool dyn, 
		const std::string& description, 
		OptionValidation< std::vector<T> >* condition,
		std::vector<T>* var) : 
    ConfigOption(name, dyn, description), 
    m_var(var),
    m_condition(condition)
  {
  }
  
  // destructor
  ~GenericOption() {}
  
  // configure the value corresponding to this option
  void resetValue(const std::string& value) 
  {
    using namespace std;
    
    // clear the array
    if (m_var->size() > 0) m_var->clear();
    
    // case of an array with one entry
    if (value.find(" ") == string::npos && 
	value.find("\t") == string::npos && 
	value.find(" \t")) {
      addEntry(value);
    }
    else { 
      size_t pos = 0;
      bool found = false;
      while (pos < value.size()) {
	size_t sep = value.find(" ", pos);
	if (sep != string::npos) {
	  addEntry(value.substr(pos, (sep-pos)));
	  pos = sep+1;
	  found = true;
	}
	
	if (sep == string::npos && (!found)){
	  sep = value.find("\t", pos);
	  if (sep != string::npos) {
	    addEntry(value.substr(pos, (sep-pos)));
	    pos = sep+1;
	    found = true;
	  }
	}	
	
	if (sep == string::npos && (!found)){
	  sep = value.find(" \t", pos);
	  if (sep != string::npos) {
	    addEntry(value.substr(pos, (sep-pos)));
	    pos = sep+1;
	    found = true;
	  }
	}
	
	// in case the separator is not found but previous 
	// entries have been detected, store the last one
	if (sep == string::npos && found) {
	  addEntry(value.substr(pos, (value.size()-pos)));
	  pos = value.size();
	}
      }
    }
  }
  
  // get the double or int data
  void getData(double* realData, int* intData) 
  {
    GetDataFun< std::vector<T> >(realData, intData, m_var);
  }
  
  // validate the option
  bool isValid() 
  {
    if (m_condition.get() != NULL) {
      return m_condition->isValid(*m_var);
    }
    return true;
  }
   
  // print value
  void print(std::ofstream& out) const 
  {
    for (unsigned i = 0; i < m_var->size(); ++i) {
      out << (*m_var)[i] << " ";
    }
    out << std::endl;
  }
  
private:
  
  // add an entry in the vector
  void addEntry(std::string entry) 
  {  
    trim(entry);	
    if (entry.size() > 0) { 
      // std::cout << " entry ###" << entry << "###" <<std::endl;
      T tmp = T();
      ConvertFromString<T>(&tmp, entry);
      m_var->push_back(tmp);
    }
  }
  
private:
  
  // variable to configure
  std::vector<T>* m_var; 
  
  // option validity condition
  SharedPtr< OptionValidation< std::vector<T> > > m_condition;
  
};

//----------------------------------------------------------------------------//
// This class implements a generic option to be configured automatically
// @author Andrea Lani

template <>
class GenericOption<bool> : public ConfigOption {
public:
  
  // constructor
  GenericOption(const std::string& name, bool dyn, 
		const std::string& description, 
		OptionValidation<bool>* condition, 
		bool* var) : 
    ConfigOption(name, dyn, description), 
    m_var(var),
    m_condition(condition)
  {
  }
  
  // destructor
  ~GenericOption() {}
  
  // get the double or int data
  void getData(double* realData, int* intData) 
  {
    *intData = (*m_var) ? 1 : 0;
  }
  
  // configure the value corresponding to this option
  void resetValue(const std::string& value) 
  {
    // convert string to value
    if (value == "true"  || value == "1") *m_var = true;
    if (value == "false" || value == "0") *m_var = false;
  }
  
  // validate the option
  bool isValid() 
  {
    if (m_condition.get() != NULL) {
      assert(m_var != NULL);
      return m_condition->isValid(*m_var);
    }
    return true;
  }
   
  // print value
  void print(std::ofstream& out) const 
  {
    const std::string flag = (*m_var) ? "true" : "false";
    out << flag << std::endl;
  }
  
private:
  
  // variable to configure
  bool* m_var; 
  
  // option validity condition
  SharedPtr< OptionValidation<bool> > m_condition;
  
};

//----------------------------------------------------------------------------//

} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
