// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoDc_hh
#define ShockFitting_CoDc_hh

//----------------------------------------------------------------------------//

#include <cmath>
#include <string>
#include "MathTools/Array2D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CoDc, whose task is to compute discontinuity
/// points by using Newthon Raphson method

class CoDc {
public:

  /// Constructor
  CoDc();

  /// Destructor
  ~CoDc();

  /// compute the new doscontinuity points values
  void callCoDc(std::vector<double>, std::vector<double>);

  /// return new shock upstream values
  std::vector<double> getnewDownValues() const {return x1; }

  /// return new downstream values
  std::vector<double> getnewUpValues() const { return x2; }

  /// return new value of the discontinuity speed
  double getnewDiscSpeed() const { return W; }

private: // helper functions

  /// assign variables used in CoDc to PhysicsData pattern
  void setPhysicsData();

  double fDc(unsigned, std::vector<double>);

private: // data

  /// gas model
  std::vector<std::string>* model;

  /// specific heat ratio
   double* gam;

  /// @param gm1 = (*gam)-1
  double gm1;

  /// vector stores downstream state
  std::vector<double>x1;

  /// vector stores upstream state
  std::vector<double>x2;

  /// Riemann invariants
  double R1, R2, S1, S2;

  /// discontinuity speed
  double W;

  /// variables used to compute the solutions
  double dyn1, dum, dum1, dum2;
  std::vector <double> dyn;
  std:: vector <double> yn;
  std:: vector <double> yn1;
  std::vector <double> b;
  Array2D <double> g;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoDc_hh

