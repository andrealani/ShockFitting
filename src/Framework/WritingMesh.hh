// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WritingMesh_hh
#define ShockFitting_WritingMesh_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"
#include "SConfig/Provider.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a WritingMesh, whose task is to write mesh and 
/// physics data computed by the code in a mesh generator format

class WritingMesh : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<WritingMesh> PROVIDER;

  /// Constructor 
  /// @param objectName the concrete class name
  WritingMesh(const std::string& objectName);

  /// Destructor
  virtual ~WritingMesh();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// write one or more output files
  virtual void write() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "WritingMesh";}

protected: // functions

  // get the name of the parent
  std::string getParentName() const {return getClassName();}

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_WritingMesh_hh
