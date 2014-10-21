// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Burgy_hh
#define ShockFitting_Burgy_hh

//--------------------------------------------------------------------------//

#include <cmath>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines Burgy, whose task is 

class Burger {
public:

  /// Constructor
  Burger();

  /// Destructor
  ~Burger();

  /// 
  double getBurg(double x)
  {
    const double UR = 5;
    double help = exp(UR/2);
    UL = UR*help/(help-1);
    double xs = 0;

    if(x<=xs) { return UL; }
    else      { return UL*(x-1); }
  }

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Burgy.hh
