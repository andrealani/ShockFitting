// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_PhysicsConsts_hh
#define ShockFitting_PhysicsConsts_hh

//--------------------------------------------------------------------------//

#include<cmath>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// Provides an a set of static functions for physics constants.

class PhysicsConsts {
public:

  /// Prandtl number 
  static double Pr() {return 0.72;}

  /// first constant of the sutherland law
  static double SuthC1() {return 1.458*pow(10,-6);}

  /// second constant of the sutherland law
  static double SuthC2() {return 110.4;}
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // Shockfitting_PhysicsConsts
