// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Prim2ParamPgDimensional_hh
#define ShockFitting_Prim2ParamPgDimensional_hh

//--------------------------------------------------------------------------//

#include "VariableTransformerSF/Prim2Param.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Prim2ParamPgDimensional, whose task is to transform
/// primitive variables in dimensional Roe parameter vector variables for a
/// perfect gas model.

class Prim2ParamPgDimensional : public Prim2Param {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Prim2ParamPgDimensional(const std::string& objectName);

  /// Destructor
  virtual ~Prim2ParamPgDimensional();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// transform variables
  virtual void transform();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Prim2ParamPgDimensional_hh