// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_PhysicsData_hh
#define ShockFitting_PhysicsData_hh

//--------------------------------------------------------------------------//

#include <map>
#include "SConfig/ConfigObject.hh"
#include "SConfig/ConfigFileReader.hh"
#include "SConfig/SharedPtr.hh"

#include "Framework/Log.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

class PhysicsData : public SConfig::ConfigObject {
public:

  /// this can be called using PhysicsData::getInstance().
  static PhysicsData& getInstance() {static PhysicsData pd; return pd;}
  
  /// Set up this object before its first use
  virtual void setup()
  {
    m_phInfo.ptr()->setup();
    m_chInfo.ptr()->setup();
    m_refInfo.ptr()->setup();
  }
  
  /// Unset up this object before its first use
  virtual void unsetup() 
  {
    m_phInfo.ptr()->unsetup();
    m_chInfo.ptr()->unsetup();
    m_refInfo.ptr()->unsetup();
  }
  
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

  /// create data with a give type and name
  template <typename ARRAYPD>
  void createData(const std::string& name)
  {
   ARRAYPD* dataArray = new ARRAYPD();
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

  /// get the physics info
  SConfig::SharedPtr<PhysicsInfo> getPhysicsInfo() const
  {
    return m_phInfo.ptr();
  }

  /// get the chemical info
  SConfig::SharedPtr<ChemicalInfo> getChemicalInfo() const 
  {
    return m_chInfo.ptr();
  }

  /// get the reference info
  SConfig::SharedPtr<ReferenceInfo> getReferenceInfo() const
  {
    return m_refInfo.ptr();
  }  

private:

  PhysicsData() : SConfig::ConfigObject("PhysicsData")
  {
    m_phInfo.name() = "PhysicsInfo";
    addOption("PhysicsInfo",&m_phInfo, "Physics info object");

    m_chInfo.name() = "ChemicalInfo";
    addOption("ChemicalInfo",&m_chInfo, "Chemical info object");

    m_refInfo.name() = "ReferenceInfo";
    addOption("ReferenceInfo",&m_refInfo, "Reference info object");
  } 

  PhysicsData(PhysicsData&) : SConfig::ConfigObject("PhysicsData") {};

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix)
  {
    LogToScreen(VERBOSE, "PhysicsData::configure() => start\n");
    SConfig::ConfigObject::configure(cmap, prefix);

    // create the physics info object
    m_phInfo.ptr().reset(new PhysicsInfo("PhysicsInfo"));
    configureDeps(cmap,m_phInfo.ptr().get());

    // create the chemical info object 
    m_chInfo.ptr().reset(new ChemicalInfo("ChemicalInfo"));
    configureDeps(cmap,m_chInfo.ptr().get());

    // create the reference info object
    m_refInfo.ptr().reset(new ReferenceInfo("ReferenceInfo"));
    configureDeps(cmap,m_refInfo.ptr().get());

    LogToScreen(VERBOSE, "PhysicsData::configure() => end\n");
  }  
  
private: // data

  /// SharePtr to PhysicsInfo object
  PAIR_TYPE(PhysicsInfo) m_phInfo;

  /// SharedPtr to ChemicalInfo object
  PAIR_TYPE(ChemicalInfo) m_chInfo;

  /// SharedPtr to ReferenceInfo object
  PAIR_TYPE(ReferenceInfo) m_refInfo;

  std::map<std::string, void*> m_mapName2ArrayPD;
  
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_PhysicsData_hh
