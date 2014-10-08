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

/// This class defines a StandardShockFitting, whose task is to 
///

class StandardShockFitting : public ShockFittingObj {
public:

  ///Constructor
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

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputFile3;
  
  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputValues1;

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputValues2;

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputValues3;

  /// command object setting bonudary node pointers
  SConfig::SharedPtr<Remeshing> m_BndryNodePtr;

  /// command object redistribute shock points
  SConfig::SharedPtr<Remeshing> m_RedistrShockPoints; 
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif //ShockFitting_StandardShockFitting_hh
