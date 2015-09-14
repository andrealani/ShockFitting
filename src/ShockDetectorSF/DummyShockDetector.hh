// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyShockDetector_hh
#define ShockFitting_DummyShockDetector_hh

//--------------------------------------------------------------------------//

#include "Framework/ShockDetector.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyShockDetector, whose task is to find the
/// location of shock discontinuities within a CFD solution.

class DummyShockDetector : public ShockDetector {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DummyShockDetector (const std::string& objectName);

  /// Destructor
  virtual ~DummyShockDetector();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command finding discontinuities 
  virtual void detect();

  /// command finding discontinuities with input/output arguments
  /// @param firstInfoVector   vector storing info used in the shock
  ///                          detector algorithm
  virtual void detect(std::vector<double>& firstInfoVector);
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_DummyShockDetector_hh
