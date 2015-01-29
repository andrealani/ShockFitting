// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MeshData_hh
#define ShockFitting_MeshData_hh

//--------------------------------------------------------------------------//

#include <map>
#include "SConfig/ConfigObject.hh"
#include "SConfig/ConfigFileReader.hh"
#include "SConfig/SharedPtr.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

class MeshData : public SConfig::ConfigObject {
public:

  /// this can be called using MeshData::getInstance().
  static MeshData& getInstance() {static MeshData md; return md;}

   /// Set up this object before its first use
  virtual void setup()
  { 
    setMeshData();
    if ( m_naddholes != 0 ) {
     caddholes->resize( 2 * m_naddholes );
     for (unsigned i=0; i<caddholes->size(); i++) {
     caddholes->at(i) = m_caddholes.at(i);}
    }
    else {caddholes->resize(1); caddholes->at(0)=m_caddholes.at(0);}
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

  template <typename ARRAYMD>
  void createData(const std::string& name)
  {
   ARRAYMD* dataArray = new ARRAYMD();
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

  /// get the number of starting step
  unsigned getnbBegin() const { return m_nbegin; }

  /// get the number of steps
  unsigned getnbSteps() const { return m_nsteps; }

  /// get the number of processors
  unsigned getnbProcessors() const { return m_nproc; }

  /// set the i-step
  void setIstep(unsigned istep) { m_istep = istep; }

  /// get the i-step
  unsigned getIstep() const { return m_istep; }

  /// get number of iterations before saving solution
  unsigned getnbIbak() const { return m_ibak; }

  /// get distance between two shock faces
  double getEPS() const { return m_eps; }

  /// get the max non dimensional distance of phantom nodes
  double getSNDMIN() const { return m_sndmin; }

  /// get the length of the shock edges
  double getDXCELL() const {return m_dxcell; }

  /// get the relax coefficient of the shock points integration
  double getSHRELAX() const { return m_shrelax; }

  /// get the bumber of additional holes
  unsigned getnbAddHoles() const { return m_naddholes; }

  /// get the Shock Fitting Version
  std::string getVersion() const { return m_version; }

  /// set the Shock Fitting Version
  void setVersion(std::string i_version) { m_version = i_version; }
 
  /// get the COOLFluiD output file
  bool withP0() { return m_coolfluidCFmesh; }

  /// get the grid informations
  bool cellsFreezed() { return m_gridInfo; }

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
    m_gridInfo = false;
    addOption("freezedWallCells",&m_gridInfo,
              "Cell close to the wall are freezed or not freezed");
    m_coolfluidCFmesh = true;
    addOption("WithP0",&m_coolfluidCFmesh,
             "Coolfluid output file");
    m_nproc = 1;
    addOption("NPROC",&m_nproc,
             "Number of processor");
    m_nbegin = 0;
    addOption("NBegin",&m_nbegin,
              "Number of starting step");
    m_nsteps = 1;
    addOption("NSteps",&m_nsteps,
              "Number of steps");
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
    caddholes = 
      MeshData::getInstance().getData < std::vector<double> > ("CADDholes");
  }

private: // data

  std::map<std::string, void*> m_mapName2ArrayMD;

private: // data
 
  /// Shock Fitting version 
  std::string m_version;

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

  /// number of strating step
  unsigned m_nbegin;

  /// number of steps
  unsigned m_nsteps;

  /// number of hole points
  unsigned m_naddholes;

  /// coolfluid output file
  bool m_coolfluidCFmesh;

  /// specifies freezed grids
  bool m_gridInfo;

  /// hole points coordinates
  std::vector <double> m_caddholes; 

  /// number of processor
  unsigned m_nproc;

  /// current executing step
  unsigned m_istep;

  /// hole points coordinates
  /// (assignable to MeshData)
  std::vector <double>* caddholes;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_MeshData_hh
