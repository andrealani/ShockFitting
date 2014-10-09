// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MeshData_hh
#define MeshData_hh

#include <map>
#include "SConfig/ConfigObject.hh"
#include "Framework/Log.hh"

class MeshData : public SConfig::ConfigObject{
public:

  /// this can be called using MeshData::getInstance().
  static MeshData& getInstance() {static MeshData md; return md;}

  /// get data with a given size, type and name
  template <typename ARRAYMD>
  ARRAYMD* getData(const std::string& name)
  {
   return (ARRAYMD*)(m_mapName2ArrayMD.find(name)->second);
  }
 
  /// create data with a given size, type and name
  template <typename ARRAYMD>
  void createData(const std::string& name, const unsigned dataSize)
  {
   ARRAYMD* dataArray = new ARRAYMD(dataSize);
   m_mapName2ArrayMD[name] = (void*)(dataArray);
  }

  /// delete data with a given size, type and name
  template <typename ARRAYMD>
   void deleteData(const std::string& name)
  {
    delete getData<ARRAYMD>(name);
  }
  
  /// get the name of the parent
  std::string getParentName() const {return std::string("MeshData");}
  
protected:
 
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix)
  {
    LogToScreen(VERBOSE, "MeshData::configure() => start\n");
    SConfig::ConfigObject::configure(cmap, prefix);
    LogToScreen(VERBOSE, "MeshData::configure() => end\n");
  }
  
private:
  MeshData() : SConfig::ConfigObject("MeshData") {}
  MeshData(MeshData&) : SConfig::ConfigObject("MeshData") {}
  
  std::map<std::string, void*> m_mapName2ArrayMD;
};

#endif //MeshData_hh
