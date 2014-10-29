// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Jcycl_hh
#define Jcycl_hh

//--------------------------------------------------------------------------//

#include <cmath>

//--------------------------------------------------------------------------//

/// This class define Jcycl, whose task is bring the value of I back 
/// into the interval [1,3] in a cyclic way.

class Jcycl {
public:

  /// Constructor
  Jcycl() {};

  /// Destructor
  ~Jcycl() {};

  /// return the value of I
  int callJcycl (unsigned I) {  int IM;
                                IM = I%3;
                                return IM+3*((1-copysign(1,(IM-1)))/2); }

};

#endif // ShockFitting_Icycl_hh
