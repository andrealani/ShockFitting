// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// //
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Area_hh
#define Area_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include <vector>

//--------------------------------------------------------------------------//

/// This class defines an Area, whose task is to compute the area bounded 
/// by the closed polygonal curve which passes through the points in 
/// the specified order with given a sequence of NB points in the plane.
/// Each simple closed curve is positively oriented (bounds positive area)
/// if and only if the points are specified in the counterclockwise order.
/// The last point of the curve is taken to be the first point specified,
/// and thus this point need not to be specified twice. However, any point
/// may be specified more than one in order to define a multiply connected
/// domain.
/// The area of triangulation may be computed by calling area with values
/// NB and NODES determined by subroutine BNODES.

class Area {
public:

  /// Constructor
  Area();

  /// Destructor
  ~Area();

  /// compute area
  void computeArea(std::vector<double>, std::vector<double>,
                   std::vector<unsigned>);

  /// return @param area
  double getArea() const { return area; }

private: // data

  /// signed area bounded by the polygonal curve defined
  /// above.
  double area;

  /// N-vectors of coordinates of points in the plane for
  /// N >= NB. NODE-I has coordinates x[I],y[I] for
  /// I=0,1,2,..N-1
  std::vector<double> x;
  std::vector<double> y;

  /// vector of node indices in the range 1 to N defining
  /// the polygonal curve
  std::vector<unsigned> nodes;

};

//--------------------------------------------------------------------------//

#endif // Area_hh
