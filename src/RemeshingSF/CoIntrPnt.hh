// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoIntrPnt_hh
#define ShockFitting_CoIntrPnt_hh

//----------------------------------------------------------------------------//

#include <vector>

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

class CoIntrPnt {
public:

  /// Constructor
  CoIntrPnt();

  /// Destructor
  ~CoIntrPnt();

  /// return x
  double getX() const;

  /// return y
  double getY() const;

  /// solve algebric system
  void callCoIntrPnt(std::vector <double>, std::vector <double>,
		     std::vector <double>, std::vector <double> );

private: // data

  /// vector of the solution variables
  std::vector <double> x;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoIntrPnt_hh
