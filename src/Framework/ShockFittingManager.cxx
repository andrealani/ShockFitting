// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ShockFittingManager.hh"
#include "Framework/Log.hh"
#include "SConfig/ConfigFileReader.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

ShockFittingManager::ShockFittingManager() : Counter(), ConfigObject("")
{
  m_ShockFittingObj.name() = "ShockFittingObj";
  addOption("ShockFittingObj",&m_ShockFittingObj.name(), "Name of the coupling tool object"); 
}
  
//--------------------------------------------------------------------------//
  
ShockFittingManager::~ShockFittingManager()
{
}
  
//--------------------------------------------------------------------------//
  
void ShockFittingManager::configure(const std::string& input, int argc, char** argv)
{   
  LogToScreen(VERBOSE, "FileProcessing::configure() root => start\n");

  // read the input file
  ConfigFileReader::OptionMap cmap;
  ConfigFileReader cfileReader;
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => before reading from " << input << "\n");
  
  // read the input file
  cfileReader.read(input, cmap);
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => after reading from " << input << "\n");
  
  const string prefix = ""; 
  configure(cmap, prefix);
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => after configure\n");
  
  // dump the options which haven't been used
  ConfigMap::getInstance().dumpWrongOptions(cmap);
  ConfigMap::getInstance().dumpOptions(); 
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => after dump options\n");
  
  // clean up ConfigFileReader::OptionMap
  ConfigFileReader::OptionMap::MapItr it;
  for (it = cmap.begin(); it != cmap.end(); ++it) {
    // manually delete the bool* which could create memory leaks
    if (it->second.second != NULL) {delete it->second.second;}
  }
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => after deallocation\n");
  
  // parse command line options
  for (int i = 0; i < argc; ++i) {
    if (string(argv[i]) == "--log") {
      int level = 0;
      istringstream iss(argv[i+1]); iss >> level;
      LogSingleton::getInstance().setLevel((LogLevel)level);
      LogToScreen(VERBOSE, "Log level is " << level << "\n");
    }  
  }
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() root => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingManager::configure(OptionMap& cmap, const string& prefix)
{
  LogToScreen(VERBOSE, "ShockFittingManager::configure() => start\n");
  
  ConfigObject::configure(cmap, prefix);
  if (ConfigFileReader::isFirstConfig()) {
    const string name = m_ShockFittingObj.name(); 
    LogToScreen(VERBOSE, "ShockFittingManager::configure() => name for coupling tools obj is " << name <<"\n");
    m_ShockFittingObj.ptr().reset(SConfig::Factory<ShockFittingObj>::getInstance().
				   getProvider(name)->create(name));
  }
  configureDeps (cmap, m_ShockFittingObj.ptr().get());
  
  LogToScreen(VERBOSE, "ShockFittingManager::configure() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingManager::setup()
{
  LogToScreen(VERBOSE, "ShockFittingManager::setup() => start\n");
  
  m_ShockFittingObj.ptr()->setup();
  
  LogToScreen(VERBOSE, "ShockFittingManager::setup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingManager::unsetup()
{
  LogToScreen(VERBOSE, "ShockFittingManager::unsetup() => start\n");
  
  m_ShockFittingObj.ptr()->unsetup();
  
  LogToScreen(VERBOSE, "ShockFittingManager::unsetup() => end\n");
}
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting
