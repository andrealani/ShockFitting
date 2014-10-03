// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
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

/// This class defines a ReadInterpolateWrite, whose task is to read a file, 
/// remesh a field and write another file with the interpolated field
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
  SConfig::SharedPtr<MeshGenerator> m_readInputFile;

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputValues1;

  /// command object reading mesh generator files
  SConfig::SharedPtr<MeshGenerator> m_readInputValues2;
  
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif //ShockFitting_StandardShockFitting_hh
