// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockFittingManager_hh
#define ShockFitting_ShockFittingManager_hh

//--------------------------------------------------------------------------//

#include "SConfig/SharedPtr.hh"
#include "SConfig/ConfigObject.hh"
#include "Framework/ShockFittingObj.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {
  
//--------------------------------------------------------------------------//

/// This class defines a ShockFittingManager, whose task is to manage coupling tools
/// 
/// @author Andrea Lani

class ShockFittingManager : public SConfig::Counter, 
			       public SConfig::ConfigObject {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  ShockFittingManager();
  
  /// Destructor
  ~ShockFittingManager();
  
  /// Set up this object before its first use
  void setup();
  
  /// Unset up this object after its last use
  void unsetup();
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  void configure(const std::string& input, int argc, char** argv);
  
  /// get the name of the parent
  std::string getParentName() const {return getName();}
  
  /// get the coupling tools object
  ShockFittingObj& getObj() {return *m_ShockFittingObj.ptr();}

protected:  
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
protected:
  
  /// coupling tools object
  PAIR_TYPE(ShockFittingObj) m_ShockFittingObj; 
  
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ShockFittingManager_hh
