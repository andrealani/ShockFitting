// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyCFDSolver_hh
#define ShockFitting_DummyCFDSolver_hh

//--------------------------------------------------------------------------//

#include "Framework/CFDSolver.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyCFDSolver, whose task is call a CFD solver.

class DummyCFDSolver : public CFDSolver {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DummyCFDSolver(const std::string& objectName);

  /// Destructor
  virtual ~DummyCFDSolver();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// Read input files
  virtual void call();
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_DummyCFDSolver_hh
