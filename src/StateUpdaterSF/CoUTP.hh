// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoUTP_hh
#define ShockFitting_CoUTP_hh

//----------------------------------------------------------------------------//

#include <cmath>
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CoUTP, whose task is to compute unknown states of the
/// triple point scheme using Newton-Raphson method.

class CoUTP {
public:

  /// Constructor
  CoUTP();

  /// Destructor
  ~CoUTP();

  /// compute unknows variables
  void callCoUTP(std::vector<double>, double, double,
                 double, double, double);

  /// return computed unknows variables
  std::vector<double> getState() const { return yn; }

  /// return exact Riemann invariant
  double getR() const { return R14; }

  /// return boolean variable
  bool getiFail() const { return ifail; }

private: // helper functions

  /// return class name
  std::string getClassName() const { return std::string("CoUTP"); }

  /// assign variables used in CoUTP to PhysicsData pattern
  void setPhysicsData();

  double fUTP(std::vector<double>, unsigned);

private: // data

  /// specific heat ratio (this value is read in PhysicsInfo)
  double* gam;

  /// matrix storing variables for the Newton-Raphon method
  Array2D <double> a;

  /// vector storing variables for the Newton-Raphon method
  std::vector <double> b;

  /// dummy variables using for the system solution
  double dum ,dumold, dum1, dum2, dyn1;
  std::vector<double> dyn;

  /// computed vector of the unknows variables
  std::vector<double> yn;

  /// work vector storing unknowns variables
  std::vector<double> yn1;

  /// Riemann invariants
  double R14, R23;

  /// 
  double DXR14, DYR14;

  /// incident shock normal speed
  double unSh;

  /// boolean variable checking convergence
  bool ifail;

  /// @param delta= (gam-1)/2
  double delta;

  /// store log infos
  FileLogManip logfile;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoUTP_hh
