// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PhysicsData_hh
#define PhysicsData_hh

//--------------------------------------------------------------------------//

#include <map>
#include <string>

class PhysicsData {
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


private: //data

 PhysicsData() {}
 PhysicsData(PhysicsData&) {}

 std::map<std::string, void*> m_mapName2ArrayPD;

};

#endif //PhysicsData_hh
