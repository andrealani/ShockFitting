// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_StandardShockFitting_hh
#define ShockFitting_StandardShockFitting_hh

//--------------------------------------------------------------------------//

#include "Framework/ShockFittingObj.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a StandardShockFitting, whose task is to read files
/// which informations about a mesh and a chemical enviroment and process
/// the data.

class StandardShockFitting : public ShockFittingObj {
public:

  ///Constructor
  /// @param objectName the concrete class name
  StandardShockFitting(const std::string& objectName);

  ///Destructor
  virtual ~StandardShockFitting();

  ///Set up this object before its first use
  virtual void setup();

  ///Unset up this object before its last use
  virtual void unsetup();

  ///Run the coupling tools
  virtual void process();

protected:

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private:

  /// Standard Shock Fitting Version
  std::string m_version;

  /// set the directory path for the results
  std::string m_resultsDir;

  /// specifies the starting file format
  bool m_startFiles;

  /// specifies if the shock fitting residuals will be computed
  bool m_computeShockFittingResidual;

  /// command object convert shock input file in sh00.dat
  SConfig::SharedPtr<Converter> m_createshockfile;

  /// command object convert triangle files 
  /// to start the Shock Fitting algorithm
  SConfig::SharedPtr<Converter> m_createTriangleFiles;

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputFile1;

  /// command object making mesh backup
  SConfig::SharedPtr<CopyMaker> m_meshBackup;
  
  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputFile2;

  /// command object fix the mesh connectivity
  SConfig::SharedPtr<MeshGenerator> m_setMeshConnect;

  /// command object setting bonudary node pointers
  SConfig::SharedPtr<Remeshing> m_bndryNodePtr;

  /// command object setting boundary faces
  SConfig::SharedPtr<Remeshing> m_bndryFacePtr;

  /// command object redistributing shock points
  SConfig::SharedPtr<Remeshing> m_redistrEqShockPoints; 

  /// command object finding phantom points
  SConfig::SharedPtr<Remeshing> m_findPhantPoints;

  /// command object changing boundary pointers
  SConfig::SharedPtr<Remeshing> m_changeBndryPoints;

  /// command object computing shock points normal
  SConfig::SharedPtr<Remeshing> m_computeNormalVector;

  /// command object remeshing shock layer
  SConfig::SharedPtr<Remeshing> m_computeShockLayer;

  /// command object fixing mesh around special points
  SConfig::SharedPtr<Remeshing> m_fixMeshSpecialPoints;

  /// command object computing shock fitting residual
  SConfig::SharedPtr<StateUpdater> m_computeSFresidual;

  /// command object writing output files
  SConfig::SharedPtr<WritingMesh> m_writeTriangleFile;

  /// command object write output files for freezed connectivity
  SConfig::SharedPtr<WritingMesh> m_writeTriangleFileFreezedConnect;

  /// command object generating new mesh
  SConfig::SharedPtr<MeshGenerator> m_callTriangle;

  /// command object caling triangle library without input file
  SConfig::SharedPtr<MeshGenerator> m_callTriangleLib;

  /// command object converting file format from Triangle to 
  /// coolfluid format (
  SConfig::SharedPtr<Converter> m_triangleToCFfmt;

  /// command object converting file format from Triangle to CFmesh
  /// when the connectivity is freezed
  SConfig::SharedPtr<Converter> m_triangleToCFmeshFreezedConnect;

  /// command object calling CFDSolver
  SConfig::SharedPtr<CFDSolver> m_COOLFluiD;

  /// command object converting file format from CFmesh to Triangle
  SConfig::SharedPtr<Converter> m_CFmeshToTriangle;

  /// command object copying Roe values
  SConfig::SharedPtr<CopyMaker> m_copyZRoe1_0;

  /// command object updating solution
  SConfig::SharedPtr<StateUpdater> m_updateSolution;

  /// command object fixing mesh around special points
  SConfig::SharedPtr<StateUpdater> m_fixSpecPoints;

  /// command object copying Roe values
  SConfig::SharedPtr<CopyMaker> m_copyZRoeSh0_1;

  /// command object moving shock points
  SConfig::SharedPtr<StateUpdater> m_moveShPoints;

  /// command object updating values in the phantom nodes
  SConfig::SharedPtr<StateUpdater> m_updatePhantPoints;

  /// command object redistributing shock points
  SConfig::SharedPtr<Remeshing> m_redistrShockPoints;

  /// command object writing back triangle node file
  SConfig::SharedPtr<WritingMesh> m_writeBackTriangleFile;

  /// command object writing shock infos
  SConfig::SharedPtr<WritingMesh> m_writeShockInfo;

  /// command object restoring mesh arrays
  SConfig::SharedPtr<CopyMaker> m_meshRestore;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_StandardShockFitting_hh
