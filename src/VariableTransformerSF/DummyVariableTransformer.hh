// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyVariableTransformer_hh
#define ShockFitting_DummyVariableTransformer_hh

//--------------------------------------------------------------------------//

#include "Framework/VariableTransformer.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyVariableTransformer, whose task is to transform a set 
/// of logically coupled variables (states or coordinates) into another set of 
/// variables. 
/// 

class DummyVariableTransformer : public VariableTransformer {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  DummyVariableTransformer(const std::string& objectName);
  
  /// Destructor
  virtual ~DummyVariableTransformer();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Transform one set of variable into another
  virtual void transform ();

  /// Transform one given set of variable into another
  virtual void transform (std::vector <double>&, std::vector <double>&,
                          std::vector <double>&);
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_DummyVariableTransformer_hh
