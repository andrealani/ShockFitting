// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeStateDps4TCneq_hh
#define ShockFitting_ComputeStateDps4TCneq_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "StateUpdaterSF/ComputeStateDps.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeStateDps4TCneq, whose task is to update
/// the solution within the grid-points located on the discontinuities
/// by enforcing Rankine-Hugoniot relations between the upstream and 
/// the downstream states for a TCneq model

class ComputeStateDps4TCneq : public ComputeStateDps {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ComputeStateDps4TCneq(const std::string& objectName);

  /// Destructor
  ~ComputeStateDps4TCneq();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// Update solution
  void update();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("ComputeStateDpsTCneq");}

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

  /// upstream value of vibrational energy
  double evu;

  /// downstream value of vibrational energy
  double evd;

  /// work vector  used to store upstream values
  std::vector<double> xu;

  /// working vector used to store downstream values
  std::vector<double> xd;

  /// working vector used to store upstream chemical species concentrations
  std::vector<double> alphau;

  /// working vector used to store downstream chemical species concentrations
  std::vector<double> alphad;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeStateDps4TCneq_hh
