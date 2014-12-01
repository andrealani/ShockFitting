// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeStateDps4Pg_hh
#define ShockFitting_ComputeStateDps4Pg_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "StateUpdaterSF/ComputeStateDps.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeStateDps4Pg, whose task is to update 
/// the solution within the grid-points located on the discontinuities
/// by enforcing Rankine-Hugoniot relations between the upstream and 
/// the downstream states for a perfect gas model

class ComputeStateDps4Pg : public ComputeStateDps {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ComputeStateDps4Pg(const std::string& objectName);

  /// Destructor
  ~ComputeStateDps4Pg();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// Update solution
  void update();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("ComputeStateDpsPg");}

  /// upload downstream status
  void recoverDownState(unsigned, unsigned);

  /// upload upstream status
  void recoverUpState(unsigned, unsigned);

  /// save old downstream status
  void saveDownState(unsigned, unsigned);

  /// compute new upstream or downstream status
  /// and store solution in Zroe arrays
  void computeDownState(unsigned, unsigned);

  void computeUpState(unsigned, unsigned);

private: // data

  /// work vector  used to store upstream values
  std::vector<double> xu;

  /// working vector used to store downstream values
  std::vector<double> xd;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeStateDps4Pg_hh
