// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Param2PrimPgDimensional_hh
#define ShockFitting_Param2PrimPgDimensional_hh

//--------------------------------------------------------------------------//

#include "VariableTransformerSF/Param2Prim.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Param2PrimPgDimensional, whose task is to transform Roe
/// parameter vector variables in primitive dimensional variables for 
/// a perfect gas model.
/// The new values are stored in new arrays to not overwrite the old mesh status.
/// These new values are pushed back at the end of the arrays of the old mesh
/// in the fortran version these new arrays are referred to index "1"
/// (ex: ZROE(1))

class Param2PrimPgDimensional : public Param2Prim {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Param2PrimPgDimensional(const std::string& objectName);

  /// Destructor
  virtual ~Param2PrimPgDimensional();

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

#endif // ShockFitting_Param2PrimPgDimensional_hh
