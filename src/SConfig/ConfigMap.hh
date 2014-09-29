// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ConfigMap_hh
#define ConfigMap_hh

#include <fstream>
#include <map>

#include "GenericOption.hh"
#include "ConfigException.hh"
#include "ConfigFileReader.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

typedef MMap<std::string, std::pair<std::string, bool*> > OptionMap; 

// This singleton class holds all options defined within this simulation
// @author Andrea Lani
class ConfigMap {
public:
  
  // get the instance of the class
  static ConfigMap& getInstance() {static ConfigMap cmap; return cmap;}
  
  // add an option corresponding to the 
  template <typename TYPE>
  void add(const std::string& name, const std::string cname1, 
	   const std::string& cname2, const std::string& desc,
	   bool isDynamic, OptionValidation<TYPE>* condition, 
	   TYPE* var, std::vector<ConfigOption*>& options)
  {    
    ConfigOption* op1 = new GenericOption<TYPE>(cname1, isDynamic, desc, condition, var);
    m_name2option[cname1] = op1; // store a pointer to the option locally
    options.push_back(op1); 
    
    // this allows to specify base class name instead fof concrete class
    // for options which are inherited
    if (cname2 != cname1) {
      ConfigOption* op2 = new GenericOption<TYPE>(cname2, isDynamic, desc, condition, var);
      m_name2option[cname2] = op2; // store a pointer to the option locally
      options.push_back(op2);  
    }
    
    // store a mapping between the name and the config name (1) 
    m_name2cname1[name] = cname1;
  }
  
  // get the data corresponding to the given string
  void getData(const std::string& name, 
	       double* realData, int* intData) throw()
  {
    if (m_name2cname1.count(name) > 0) {
      const std::string cname = m_name2cname1.find(name)->second;
      m_name2option.find(cname)->second->getData(realData, intData);
    }
    else throw ConfigException(std::string("ConfigMap::getData() => Option < " + 
					   name + " > is not available!"));
  }
  
  // dump options description to file
  void dumpOptions() const 
  {
    using namespace std;
    
    int rank = 0;
#ifdef USING_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    if (rank == 0) {
      ofstream fout("options.dat");
      map<string, ConfigOption*>::const_iterator it;
      for (it = m_name2option.begin(); it != m_name2option.end(); ++it) {
	fout << it->first << " : " << it->second->description() << " => ";
	it->second->print(fout);
      }
    }
  }
  
  // dump all wrongly declared options
  void dumpWrongOptions(OptionMap& omap) throw()
  {
    using namespace std;
    
    int rank = 0;
#ifdef USING_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    if (rank == 0 && ConfigFileReader::isFirstConfig()) {
      int count = 0;
      OptionMap::MapItr it;
      for (it = omap.begin(); it != omap.end(); ++it) {
	if (!(*it->second.second)) {
	  cout << "Key < " <<  it->first << " > doesn't exist\n";
	  count++;
	}
      }
      
      if (count > 0) throw ConfigException
	(string("ConfigMap::dumpWrongOptions() => Fix the errors listed above"));
    }
  }
  
  // validate all the options
  void validateOptions() throw()
  {
    using namespace std;
    
    int count = 0;
    std::map<std::string, ConfigOption*>::const_iterator it;
    for (it = m_name2option.begin(); it != m_name2option.end(); ++it) {
      if (!(it->second->isValid())) {
	cout << "Condition on Key < " << it->first 
	     << " > and related parameters (if any) is not satisfied\n";
	count++;
      }
    }
    
    if (count > 0) throw ConfigException
      (string("ConfigMap::validateOptions() => Fix the errors listed above"));
  }
  
private:
  
  // make this class not instantiable
  ConfigMap() {}
  
  ~ConfigMap() 
  { 
    using namespace std;
    
    // deallocate the pointers to ConfigOption
    map<string, ConfigOption*>::const_iterator it;
    for (it = m_name2option.begin(); it != m_name2option.end(); ++it) {
      delete it->second;
    }
  }
  
  // make this class not copyable
  ConfigMap(const ConfigMap& c);
  const ConfigMap& operator=(const ConfigMap& c);
  
private:
  
  // map with [key = option name, value = ConfigOption*]
  std::map<std::string, ConfigOption*> m_name2option;
  
  // map with [key = name, value = config name (1)]
  std::map<std::string, std::string> m_name2cname1;
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
