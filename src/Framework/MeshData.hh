// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

//#ifndef ShockFitting_MeshData_hh
//#define ShockFitting_MeshData_hh

#ifndef MeshData_hh
#define MeshData_hh

//--------------------------------------------------------------------------//

#include <map>
#include "SConfig/ConfigObject.hh"
#include "SConfig/ConfigFileReader.hh"
#include "SConfig/SharedPtr.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

//namespace Shockfitting {

//--------------------------------------------------------------------------//

class MeshData : public SConfig::ConfigObject {
public:

  /// this can be called using MeshData::getInstance().
  static MeshData& getInstance() {static MeshData md; return md;}

   /// Set up this object before its first use
  virtual void setup()
  { 
    setMeshData();
    *eps = m_eps;
    *sndmin = m_sndmin;
    *dxcell = m_dxcell;
    *shrelax = m_shrelax;
    *ibak = m_ibak;
    *naddholes = m_naddholes;
    if ( (*naddholes) != 0 ) {
     caddholes->resize( 2 * (*naddholes) );
     for (unsigned i=0; i<caddholes->size(); i++) {
     caddholes->at(i) = m_caddholes.at(i);}
    }
    else {caddholes->resize(1); caddholes->at(0)=m_caddholes.at(0);}
    *nproc = m_nproc;

  }

  /// Unset up this object before its first use
  virtual void unsetup(){}

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

private:

  MeshData() : SConfig::ConfigObject("MeshData") 
  {
    m_eps = 0;
    addOption("EPS",&m_eps,
            "Distance between two shock faces");
    m_sndmin = 0;
    addOption("SNDMIN",&m_sndmin,
            "Max non dimensional distance of phantom nodes");
    m_dxcell = 0;
    addOption("DXCELL",&m_dxcell,
            "Length of the shock edges");
    m_shrelax = 0;
    addOption("SHRELAX",&m_shrelax,
            "Relax Coefficient of shock points integration");
    m_ibak = 1;
    addOption("IBAK",&m_ibak,
            "Number of iterations before saving solution",true);
    m_naddholes = 0;
    addOption("Naddholes",&m_naddholes,
            "Number of hole points");
    m_caddholes = std::vector <double>();
    addOption("CADDholes",&m_caddholes,
            "Holes points coordinates");
    m_nproc = 1;
    addOption("NPROC",&m_nproc,
            "Number of processor");
   }
  
   MeshData(MeshData&) : SConfig::ConfigObject("MeshData") {}

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix)
  {
    LogToScreen(VERBOSE, "MeshData::configure() => start\n");
    SConfig::ConfigObject::configure(cmap, prefix);
    LogToScreen(VERBOSE, "MeshData::configure() => end\n");
  }

private: // helper functions

  void setMeshData()
  {
    eps = MeshData::getInstance().getData <double> ("EPS");
    sndmin = MeshData::getInstance().getData <double> ("SNDMIN");
    dxcell = MeshData::getInstance().getData <double> ("DXCELL");
    shrelax = MeshData::getInstance().getData <double> ("SHRELAX");
    ibak = MeshData::getInstance().getData <unsigned> ("IBAK");
    naddholes = MeshData::getInstance().getData <unsigned> ("Naddholes");
    caddholes = MeshData::getInstance().getData < std::vector<double> > ("CADDholes");
    nproc = MeshData::getInstance().getData <unsigned> ("NPROC");
  }

private: // data

  std::map<std::string, void*> m_mapName2ArrayMD;

private: // data (read from input file)

  /// distance between two shock faces
  double m_eps;

  /// max non dimensional distance of phantom nodes
  double m_sndmin;

  /// length of the shock edges
  double m_dxcell;

  /// relax coefficient for shock points integration
  double m_shrelax;

  /// number of iterations before saving solution
  unsigned m_ibak;

  /// number of hole points
  unsigned m_naddholes;

  /// hole points coordinates
  std::vector <double> m_caddholes; 

  /// number of processor
  unsigned m_nproc;

  /// distance between two shock faces
  /// (assignable to MeshData)
  double* eps;

  /// max non dimensional distance of phantom nodes
  /// (assignable to MeshData)
  double* sndmin;

  /// length of the shock edges
  /// (assignable to MeshData)
  double* dxcell;

  /// relax coefficient for shock points integration
  /// (assignable to MeshData)
  double* shrelax;

  /// number of iterations before saving solution
  /// (assignable to MeshData)
  unsigned* ibak;

  /// number of hole points
  /// (assignable to MeshData)
  unsigned* naddholes;

  /// hole points coordinates
  /// (assignable to MeshData)
  std::vector <double>* caddholes;

  /// number of processor
  /// (assignable to MeshData)
  unsigned* nproc;

};

//--------------------------------------------------------------------------//

//} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_MeshData_hh
