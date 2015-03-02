// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoShock_hh
#define ShockFitting_CoShock_hh

//----------------------------------------------------------------------------//

#include <cmath>
#include <string>
#include <vector>
#include "MathTools/Array2D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CoShock, whose task is to compute shock points
/// by using Newton-Raphson method.

class CoShock {
public:

  /// Constructor
  CoShock();

  /// Destructor
  ~CoShock();

  /// compute new shock points values
  void callCoShock(std::vector<double>, std::vector<double>, double);

  /// return new shock points values
  std::vector<double> getnewDownValues() const { return x1;}

  /// return new value of the discontinuity speed
  double getnewDiscSpeed() const { return W; }

private: // helper functions

  /// assign variables used in CoShock to PhysicsData pattern
  void setPhysicsData();

  double fShock(unsigned, std::vector<double>);

private: // data

  /// heat specific ratio
  double* gam;

  /// vector stores downstream status
  std::vector<double>x1;

  /// vector stores upstream status
  std::vector<double>x2;

  /// Riemann invariant
  double R2;

  /// discontinuity speed
  double W;

  /// variables used to compute the solutions
  double dyn1, dum, dum1, dum2;
  std::vector <double> dyn;
  std::vector <double> yn;
  std::vector <double> yn1;
  std::vector <double> b;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoShock_hh
