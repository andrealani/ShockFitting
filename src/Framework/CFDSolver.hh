// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CFDSolver_hh
#define ShockFitting_CFDSolver_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines CFDSolver, whose task is to call a CFDSolver

class CFDSolver : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<CFDSolver> PROVIDER;

  /// Constructor 
  /// @param objectName the concrete class name
  CFDSolver(const std::string& objectName);

  /// Destructor
  ~CFDSolver();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// call the CFD solver
  virtual void call() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "CFDSolver";}

protected: // data

  /// set the option to change some CFD input file values at run time
  bool m_alterCFDinputfile;

  /// specifies the values to be changed
  std::vector<std::string> m_alterWhichValues; 
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_CFDSolver_hh
