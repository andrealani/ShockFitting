// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FieldInterpolator_hh
#define ShockFitting_FieldInterpolator_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"
#include "Framework/Connectivity.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {
  class Field;
  
//--------------------------------------------------------------------------//

/// This class defines a FieldInterpolator, whose task is to transfer a 
/// discretized field from a grid into another, in space and/or time.
/// 
/// @author Andrea Lani

class FieldInterpolator : public BaseShockFitting {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<FieldInterpolator> PROVIDER;
    
  /// Constructor 
  /// @param objectName the concrete class name
  FieldInterpolator(const std::string& objectName);
  
  /// Destructor
  virtual ~FieldInterpolator();
   
  /// Set up this object before its first use
  virtual void setup() = 0;
  
  /// Unset up this object after its last use
  virtual void unsetup() = 0;
  
  /// Interpolate from one field and data structure into another
  /// @param inField    input field for the solution to be interpolated
  /// @param outField   output field resulting from the interpolation inputs
  virtual void interpolate(Field* inField, Field* outField) = 0;
  
  /// Gets the Class name
  static std::string getClassName() {return "FieldInterpolator";}
  
protected:
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
protected:
  
  /// IDs for x,y,z in the input field
  std::vector<unsigned> m_inXyzIDs;
  
  /// IDs for x,y,z in the output field
  std::vector<unsigned> m_outXyzIDs;
  
};

//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_FieldInterpolator_hh
