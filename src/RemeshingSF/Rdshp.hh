// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Rdshp_hh
#define ShockFitting_Rdshp_hh

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines Rdshp, whose task is compute the distance between
/// a cell vertex (xc, yc) and a shock segment denoted by the two
/// two shock points (xs1, ys1, xs2, ys2)

class Rdshp {
public:

  /// Constructor
  Rdshp ();

  /// Destructor
  ~Rdshp();

  /// return  distance between a cell vertex and shock segment
  double getRdshp (double, double, double,
                   double, double, double);

private: // helper function

  /// assign global variable SNDMIN to MeshData data
  void setSNDMIN();

  /// assign external values to private data
  void resetValues(double XC, double YC, double XS1,
                   double YS1, double XS2, double YS2);

  /// compute distance between a cell vertex and shock segment
  void callRdshp ();

private: // data

  /// max normalized distance of phantom point
  double* SNDmin;

  /// coordinates of cell vertex
  double xc, yc;

  /// shock points which denote straight line
  double xs1, ys1;

  /// shock points which denote straight line
  double xs2, ys2;

  /// distance
  double rdshp;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Rdshp_hh
