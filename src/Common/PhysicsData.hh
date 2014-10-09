// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PhysicsData_hh
#define PhysicsData_hh

//--------------------------------------------------------------------------//

#include <map>
#include "SConfig/ConfigObject.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

class PhysicsData : public SConfig::ConfigObject {
public:
  /// this can be called using MeshData::getInstance().
  static PhysicsData& getInstance() {static PhysicsData pd; return pd;}

 /// get data with a given size, type and name
 template <typename ARRAYPD>
 ARRAYPD* getData(const std::string& name)
 {
   return (ARRAYPD*)(m_mapName2ArrayPD.find(name)->second);
 }

 /// create data with a given size, type and name
 template <typename ARRAYPD>
 void createData(const std::string& name, const unsigned dataSize)
 {
   ARRAYPD* dataArray = new ARRAYPD(dataSize);
   m_mapName2ArrayPD[name] = (void*)(dataArray);
 }

 /// delete data with a given size, type and name
 template <typename ARRAYPD>
 void deleteData(const std::string& name)
 {
   delete getData<ARRAYPD>(name);
 }
  
  /// get the name of the parent
  std::string getParentName() const {return std::string("PhysicsData");}

protected:
 
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix)
  {
    LogToScreen(VERBOSE, "PhysicsData::configure() => start\n");
    SConfig::ConfigObject::configure(cmap, prefix);
    LogToScreen(VERBOSE, "PhysicsData::configure() => end\n");
  }
  
private: //data

 PhysicsData() : SConfig::ConfigObject("PhysicsData") {}
  PhysicsData(PhysicsData&) : SConfig::ConfigObject("PhysicsData") {}
  
 std::map<std::string, void*> m_mapName2ArrayPD;
  
};

#endif //PhysicsData_hh
