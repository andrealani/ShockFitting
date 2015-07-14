// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BaseShockFitting_hh
#define ShockFitting_BaseShockFitting_hh

//--------------------------------------------------------------------------//

#include "SConfig/ConfigObject.hh"
#include "SConfig/Provider.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a BaseShockFitting, whose task is to convert units for 
/// individual physical quantities from one system to another
/// 
/// @author Andrea Lani

class BaseShockFitting : public SConfig::Counter, 
			 public SConfig::ConfigObject {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<BaseShockFitting> PROVIDER;
  
  /// Constructor 
  /// @param objectName the concrete class name
  BaseShockFitting(const std::string& objectName);
  
  /// Destructor
  virtual ~BaseShockFitting();
  
  /// Set up this object before its first use
  virtual void setup() = 0;
  
  /// Unset up this object after its last use
  virtual void unsetup() = 0;
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
  /// Gets the Class name
  static std::string getClassName() {return "BaseShockFitting";}

protected: // functions
  
  // get the name of the parent
  std::string getParentName() const {return getClassName();}

protected: // data

  /// flag input variables which are active
  std::vector<bool> m_flagInputVars;
    
  /// flag output variables which are active
  std::vector<bool> m_flagOutputVars;
         
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_BaseShockFitting_hh
