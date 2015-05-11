// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockFittingObj_hh
#define ShockFitting_ShockFittingObj_hh

//--------------------------------------------------------------------------//

#include "SConfig/ConfigObject.hh"
#include "SConfig/Factory.hh"

#include "Framework/Converter.hh"
#include "Framework/CopyMaker.hh"
#include "Framework/CFDSolver.hh"
#include "Framework/FieldInterpolator.hh"
#include "Framework/FileProcessing.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Remeshing.hh"
#include "Framework/StateUpdater.hh"
#include "Framework/VariableTransformer.hh"
#include "Framework/WritingMesh.hh"
#include "RemeshingSF/CoNorm.hh"
#include "StateUpdaterSF/ComputeStateDps.hh"
#include "StateUpdaterSF/MoveDps.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//--------------------------------------------------------------------------//

namespace ShockFitting {
  
  class BoundaryConnectivity;
  class Field_SF;
  
//--------------------------------------------------------------------------//

/// This class defines a ShockFittingObj, whose task is to provide an 
/// interface to the individual tools.
/// 
/// @author Andrea Lani

class ShockFittingObj : public SConfig::Counter, 
			public SConfig::ConfigObject {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<ShockFittingObj> PROVIDER;
  
  /// Constructor 
  /// @param objectName the concrete class name
  ShockFittingObj(const std::string& objectName);
  
  /// Destructor
  virtual ~ShockFittingObj();
  
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the coupling tools
  virtual void process();
  
  /// Run the coupling tools feeding the mesh data and getting them back modified 
  virtual void processField(const BoundaryConnectivity *const inBndConn, 
			    const Field *const inState,
			    const Field *const inNode,
			    BoundaryConnectivity* outBndConn, 
			    Field* outState,
			    Field* outNode);
  
  /// get the name of the parent
  virtual std::string getParentName() const {return getName();}

  /// get the field interpolators list
  std::vector<PAIR_TYPE(FieldInterpolator)>& getFieldInterpolatorList()
  {
    return m_fInterpolator;
  }
  
  /// get the file processing list
  std::vector<PAIR_TYPE(FileProcessing)>& getFileProcessingList()
  {
    return m_fProcessing;
  }

  /// get the meshgenerator list
  std::vector<PAIR_TYPE(MeshGenerator)>& getMeshGeneratorList()
  {
    return m_mGenerator;
  }

  /// get the remeshing list
  std::vector<PAIR_TYPE(Remeshing)>& getRemeshingList()
  {
    return m_fRemeshing;
  }

  /// get the writingmesh list
  std::vector<PAIR_TYPE(WritingMesh)>& getWritingMeshList()
  {
    return m_wMesh;
  }

  /// get the converters list
  std::vector<PAIR_TYPE(Converter)>& getConverterList()
  {
    return m_fConverter;
  }

  /// get the copy maker list
  std::vector<PAIR_TYPE(CopyMaker)>& getCopyMakerList()
  {
    return m_cMaker;
  }

  /// get the state updater list
  std::vector<PAIR_TYPE(StateUpdater)>& getStateUpdaterList()
  {
    return m_sUpdater;
  }

protected:
   
  /// create a list of coupling objects
  template <typename T>
  void createList(std::vector<PAIR_TYPE(T)>& obj)
  {
    for (unsigned i = 0; i < obj.size(); ++i) {
      const std::string name = obj[i].name();
      obj[i].ptr().reset(SConfig::Factory<T>::getInstance().
			 getProvider(name)->create(name));
    }
  }
   
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
    
protected:
  
  /// array of field interpolators
  std::vector<PAIR_TYPE(FieldInterpolator)> m_fInterpolator;
  
  /// array of file processing
  std::vector<PAIR_TYPE(FileProcessing)> m_fProcessing;

  /// array of mesh generator readers
  std::vector<PAIR_TYPE(MeshGenerator)> m_mGenerator ;

  /// array of field remeshing
  std::vector<PAIR_TYPE(Remeshing)> m_fRemeshing;

  /// object normal vector computing
  PAIR_TYPE(CoNorm) m_cNormalVector;

  /// array of data writing
  std::vector<PAIR_TYPE(WritingMesh)> m_wMesh;

  /// array of converter objects
  std::vector<PAIR_TYPE(Converter)> m_fConverter;

  /// array of copy maker objects
  std::vector<PAIR_TYPE(CopyMaker)> m_cMaker;

  /// array of state updater objects
  std::vector<PAIR_TYPE(StateUpdater)> m_sUpdater;

  /// object calling CFDSolver
  PAIR_TYPE(CFDSolver) m_CFDSolver;

  /// object updating solution
  PAIR_TYPE(ComputeStateDps) m_cState;

  /// object moving shock points
  PAIR_TYPE(MoveDps) m_moveDps;

private: // helper functions

  /// create MeshData variables
  void createMeshData();

  /// create PhysicsData variables
  void createPhysicsData();

  /// delete MeshData variables
  void deleteMeshData();

  /// delete PhysicsData variables
  void deletePhysicsData();

};  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ShockFittingObj_hh
