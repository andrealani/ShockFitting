// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Ishel_hh
#define Ishel_hh

//--------------------------------------------------------------------------//

#include <vector>

//--------------------------------------------------------------------------//

/// This class has members ishel1 and ishel2.
/// If both members return 0 then shock segment denoted by the
/// two shock points (xs1, ys1, xs2, ys2) crosses  the cell with
/// the vertices xc1, yc1, xc2, yc2 ,xc3, yc3.

class Ishel {
public:

  /// Default constructor
  Ishel ();

  /// Constructor with class params assigned
  Ishel (std::vector <double>, std::vector <double>,
         double, double, double, double);

  /// Destructor
  ~Ishel();

  /// This member returns 0 if the cell is crossed by the straight line
  /// passing for the two shock points. 
  /// Otherwise the function ishel1 returns 1.
  /// In particular it evaluates the sign of results obtained
  /// replacing the vertex coordinates in the equation of the straight
  /// line denoted by the two shock points.
  /// If all results have the same sign then the straight line does
  /// not cross the triangle.
  int Ishel1 ();

  /// This member is applied only to the cells where the member Ishel1
  /// returns 0;
  /// Ishel2 returns 0 if ?almost? one cell segment has the
  /// intersection point with the shock straight line enclosed
  /// between the two shock points.
  int Ishel2 ();

private: // data

  /// x-coordinates of cell vertex
  std::vector <double> xc;

  /// y-coordinates of cell vertex
  std::vector <double> yc;

  /// shock points which denote straight line
  double xs1, ys1;

  /// shock points which denote straight line
  double xs2, ys2;

};

#endif // Ishel_hh

