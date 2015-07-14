// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iostream>

#include "StringManip.hh"
#include "ConfigFileReader.hh"
#include "ConfigException.hh"

//----------------------------------------------------------------------------//

using namespace std;

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// flag telling if this is the first configuration
bool ConfigFileReader::m_firstConfig = true;

//----------------------------------------------------------------------------//

// read the configuration file
void ConfigFileReader::read(const string& filename, OptionMap& cmap)
{  
  m_cfile = filename;
  ifstream file(m_cfile.c_str());
  
  string line = "";
  string key = "";
  string value = "";
  bool append = false;
  
  // we'll use the set for checking for duplicated keys
  set<string> keyList;
  ofstream fout("config.log");
  
  try {
    int countBrace = 0;
    string prefix = "";
    vector<string> vprefix;
    
    bool commentBlock = false;
    while (getline(file,line)) {
      size_t commMark = line.find('$');
      if (commMark == 0) {
	commentBlock = !commentBlock;
	continue;
      }      
      
      if (!commentBlock) {
	// ignore lines starting with '#'
	size_t posComment = line.find('#');
	// proceed only if the character # is not the first character
	// or if # is not present at all
	if (posComment != 0 || posComment == string::npos) {
	  if (posComment != string::npos) {
	    // if # is present, trim out what follows that character
	    const string subs = line.substr(0,posComment);
	    line = subs;
	  }
	  
	  // check if a curly brace is present in the line
	  bool hasBrace = false;
	  string s(line); trim(s);
	  if (s.size() != 0) {
	    trim(line);
	    bool continuation = false;
	    // look for "{" for grouping options
	    const size_t lbrace = line.find('{');
	    string prefixTmp = "";
	    if (lbrace != string::npos) {
	      prefixTmp = line.substr(0,lbrace);
	      trim(prefixTmp);
	      vprefix.push_back(prefixTmp);
	      
	      prefix = "";
	      for (size_t is = 0; is < vprefix.size(); ++is) {
		prefix += vprefix[is];
	      }
	      
	      countBrace++;
	      // here you must have value = key or value {
	      hasBrace = true;
	      // if "{" is found, stop parsing, go to the next line
	    }
	    else {
	      // cout << "{ NOT found" << endl;
	      // if "{" is not found, continue parsing
	      // look for "}" for grouping options
	      const size_t rbrace = line.find('}');
	      if (rbrace != string::npos) {
		// cout << "} found" << endl;
		countBrace--;
		// remove the last key in the prefix 
		vprefix.pop_back();
		if (vprefix.size() == 0) {
		  hasBrace = false;
		  prefix = "";
		}
	      }
	      else {
		// if "}" is not found, continue parsing
		// case in which a format "A = B" (or "A = B \") is expected 
		// cout << "} NOT found" << endl;
		if ((!append) && (!hasBrace)) {
		  // look for the continuation character
		  // if present, remove it
		  const size_t posc = line.find('\\');
		  if (posc != string::npos) {
		    line.replace(posc, 1, " ");
		    trim(line);
		    append = true;
		    continuation = true;
		  }
		  
		  // separate the line between key and value:
		  // key   in the form "name1.name2.name3" 
		  // value is all the rest after the =
		  size_t start = 0;
		  size_t sep = line.find('=', start); //change
		  while(sep != string::npos) {
		    // key
		    string keyTmp = line.substr(start, sep-start);
		    trim(keyTmp);
		    key = prefix + keyTmp;
		    trim(key);
		    
		    // if there is ";" mark that as a starting position, else 
		    // start is assigned to string::npos
		    const size_t end = line.find(';', start);
		    value = (end != string::npos) ?
		      line.substr(sep+1, end-(sep+1)) : line.substr(sep+1);
		    trim(value);  
		    
		    if (!continuation) {
		      checkKeys(key, keyList); 
		      cmap.insert(key, pair<string,bool*>(value,new bool(false)));
		      fout << key << " <=> " << value << endl;
		      append = false; 
		      continuation = false;
		    }
		    
		    if (end != string::npos) {
		      start = end+1;
		      sep = line.find('=', start);
		    }
		    else {
		      start = sep = string::npos;
		    }
		  }
		}
		else if (append) {
		  // case in which I'm appending (the previous line ended wth \)
		  const size_t posc = line.find('\\');
		  if (posc != string::npos) {
		    append = true;
		    // remove the continuation character
		    line.replace(posc, 1, " "); 
		    trim(line);
		    value += " " + line;
		  }
		  else {
		    value += " " + line;
		    trim(value);  
		    checkKeys(key, keyList);
		    cmap.insert(key, pair<string,bool*>(value,new bool(false)));
		    append = false;  
		    
		    printError(fout, string(key + " <=> " + value));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    if (countBrace != 0) {
      throw ConfigException(string("number of \"{\"s doesn't match number of \"}\"s")); 
    }
  }
  catch(exception& e) {
    printError(cout, e.what());
  }
  
  file.close();
  cmap.sort();
  
  static int count = 0;
  if (count > 0) m_firstConfig = false;
  count++;
}

//----------------------------------------------------------------------------//

// check for duplicated keys
void ConfigFileReader::checkKeys(const std::string& key, 
				 std::set<std::string>& keyList) throw()
{
  // cout << "key = " << key << " is present " <<  keyList.count(key) << "times \n";
  
  if (keyList.count(key) == 0) {
    keyList.insert(key);
  }
  else {	
    throw ConfigException(string("Key < " + key  + " > is DUPLICATED"));
  }
  
  // take the option name, without any prefix
  const size_t pos = key.find_last_of(".");
  const string param = key.substr(pos+1, key.size() - (pos+1));
  if (keyList.count(param) == 0) {
    keyList.insert(param);
  }
  else { 
    printError(cout, string("WARNING: Parameter < " + param  + 
			    " > is DUPLICATED: is this what you really expect?"));
  }
}

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//
