// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyFieldInterpolator_hh
#define ShockFitting_DummyFieldInterpolator_hh

//--------------------------------------------------------------------------//

#include "Framework/FieldInterpolator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyFieldInterpolator, whose task is to transfer a 
/// discretized field from a grid into another, in space and/or time.
/// 
/// @author Andrea Lani

class DummyFieldInterpolator : public FieldInterpolator {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  DummyFieldInterpolator(const std::string& objectName);
  
  /// Destructor
  virtual ~DummyFieldInterpolator();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Interpolate from one field and data structure into another
  /// @param inField    input field for the solution to be interpolated
  /// @param outField   output field resulting from the interpolation inputs
  virtual void interpolate(Field* inField, Field* outField);
  
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_DummyFieldInterpolator_hh
