// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShpDpndnc_hh
#define ShockFitting_ShpDpndnc_hh

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class evaluates shock points dependecies

class ShpDpndnc {
public:

  /// Constructor
  ShpDpndnc(double, double, double,
            double, double, double,
            double, double, double);

  /// Destructor
  ~ShpDpndnc();

  /// return shp_dependence
  unsigned getDependence();

private: // helper function

  /// compute shp_dependence
  void callShpDpndnc();

private: // data

  /// shock point dependence
  unsigned shp_dependence;

  /// coordinates of i-shockpoint
  double x,y;

  /// coordinates if (i+1)-shockpoint
  double xi,yi;

  /// shock point speed components
  double ush, vsh;

  /// (i+1)-shockpoint status
  double ui,vi,ai;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ShpDpndnc_hh
