// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MeshData_hh
#define MeshData_hh

#include <map>
#include <string>

class MeshData {
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
 
  private:
   MeshData() {}
   MeshData(MeshData&) {}
  
 
   std::map<std::string, void*> m_mapName2ArrayMD;
};

#endif //MeshData_hh
