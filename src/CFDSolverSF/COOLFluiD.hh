// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_COOLFluiD_hh
#define ShockFitting_COOLFluiD_hh

//--------------------------------------------------------------------------//

#include <sstream>
#include "Framework/CFDSolver.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a COOLFluiD, whose task is to call COOLFluiD software

class COOLFluiD : public CFDSolver {
public:

  /// Constructor
  COOLFluiD(const std::string& objectName);

  /// Destructor
  virtual ~COOLFluiD();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command variables transformation
  virtual void call();

private: // data

  /// string for system command
  std::stringstream command;
};

//--------------------------------------------------------------------------//

}

//--------------------------------------------------------------------------//

#endif // ShockFitting_COOLFluiD_hh


