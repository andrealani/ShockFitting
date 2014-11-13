// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MoveDps4TCneq_hh
#define ShockFitting_MoveDps4TCneq_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "StateUpdaterSF/MoveDps.hh"
#include "Framework/FileLogManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

/// This class defines a MoveDps4TCneq, whose task is to compute the new 
/// position of the shock for a TCneq model

class MoveDps4TCneq : public MoveDps {
public:

  /// Constructor
  /// @param objectName the concrete class name
  MoveDps4TCneq(const std::string& objectName);

  /// Destructor
  ~MoveDps4TCneq();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// move the shock points
  void update();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("MoveDpsTCneq");}

private: // data

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_MoveDps4TCneq_hh
