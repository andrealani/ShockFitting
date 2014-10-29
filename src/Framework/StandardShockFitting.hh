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

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputFile1;
  
  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputFile2;

  /// command object setting bonudary node pointers
  SConfig::SharedPtr<Remeshing> m_bndryNodePtr;

  /// command object redistributing shock points
  SConfig::SharedPtr<Remeshing> m_redistrShockPoints; 

  /// command object finding phantom points
  SConfig::SharedPtr<Remeshing> m_findPhantPoints;

  /// command object changing boundary pointers
  SConfig::SharedPtr<Remeshing> m_changeBndryPoints;

  /// command object computing shock points normal
  SConfig::SharedPtr<CoNorm> m_computeNormalVector;

  /// command object remeshing shock layer
  SConfig::SharedPtr<Remeshing> m_computeShockLayer;

  /// command object fixing mesh around special points
  SConfig::SharedPtr<Remeshing> m_fixMeshSpecialPoints;

  /// command object writing output files
  SConfig::SharedPtr<WritingMesh> m_writeTriangleFile;

  /// command object generating new mesh
  SConfig::SharedPtr<MeshGenerator> m_callTriangle;

  /// command object converting file format from Triangle to CFmesh
  SConfig::SharedPtr<Converter> m_triangleToCFmesh;

  /// command object converting file format from CFmesh to Triangle
  SConfig::SharedPtr<Converter> m_CFmeshToTriangle;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif //ShockFitting_StandardShockFitting_hh
