// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

//--------------------------------------------------------------------------//

#ifndef ShockFitting_StateUpdater_hh
#define ShockFitting_StateUpdater_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"
#include "SConfig/Provider.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a StateUpdater, whose task is update the solution
/// within the grid-points located on the discontinuities by enforcing
/// Rankine-Hugoniot relations between the upstream and downstream
/// states.


class StateUpdater : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<StateUpdater> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  StateUpdater(const std::string& objectName);

  /// Destructor
  virtual ~StateUpdater();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// update status
  virtual void update() = 0;

  /// Gets the class name
  static std::string getClassName() {return "StateUpdater";}

protected: // functions

  // get the name of the parent
  std::string getParentName() const {return getClassName();}
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // Shockfitting_StateUpdater_hh
