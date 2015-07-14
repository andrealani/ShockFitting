// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CopyMaker_hh
#define ShockFitting_CopyMaker_hh

//--------------------------------------------------------------------------//

#include "SConfig/ConfigObject.hh"
#include "SConfig/Provider.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CopyMaker, whose task is to make copy
/// of the mesh variables and arrays.

class CopyMaker : public SConfig::Counter,
                  public SConfig::ConfigObject { 
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<CopyMaker> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  CopyMaker(const std::string& objectName);

  /// Destructor
  virtual ~CopyMaker();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// Copy arrays
  virtual void copy() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "CopyMaker";}
protected: //functions

  // get the name of the parent
  std::string getParentName() const {return getClassName();}
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_CopyMaker_hh
