// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyStateUpdater_hh
#define ShockFitting_DummyStateUpdater_hh

//--------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyStateUpdater, whose task is to enforce
/// Rankine-Hugoniot relations between the upstream and the downstream
/// state.

class DummyStateUpdater : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  DummyStateUpdater(const std::string& objectName);

  /// Destructor
  virtual ~DummyStateUpdater();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command solution update
  virtual void update();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting_hh

#endif // ShockFitting_DummyStateUpdater_hh
