// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ConfigFileReader_hh
#define ConfigFileReader_hh

#include <string>
#include <set>

#include "MMap.hh"

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements a parser for the file format 
// corresponding to the self-configuration capability
// @author Andrea Lani

class ConfigFileReader {
public:
  
  typedef MMap<std::string, std::pair<std::string, bool*> > OptionMap; 
  
  // constructor
  ConfigFileReader() : m_cfile("input.dat") {}
  
  // destructor
  ~ConfigFileReader() {}
  
  // read the configuration file
  void read(const std::string& filename, OptionMap& cmap);
  
  // static function to tell if it's the first reading
  static bool isFirstConfig() {return m_firstConfig;}
  
private:
  
  // check for duplicated keys
  void checkKeys(const std::string& key, std::set<std::string>& keyList) throw();
  
private:
  
  // configuration file name
  std::string m_cfile;
  
  // flag telling if this is the first configuration
  static bool m_firstConfig;
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
