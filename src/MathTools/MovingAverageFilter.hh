// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MovingAverageFilter_hh
#define MovingAverageFilter_hh

//--------------------------------------------------------------------------//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdio.h>
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

/// This class defines a MovingAverageFilter, whose task is to
/// smooth a curve profile by averaging the Y-components using the
/// ones of the closest points

/// @author Valentina De Amicis

class MovingAverageFilter {
public:

  /// Constructor
  MovingAverageFilter();

  /// Destructor
  ~MovingAverageFilter();

  /// smooth the curves
  /// @param Y                   y-coordinates of the curve
  /// @param nbSmoothingPoints   number of points used for the
  ///                            smoothing (that is the averaging)
  void curveForwardSmoothing(std::vector<double> Y,
                             int nbSmoothingPoints);

  /// return the new Y-component
  std::vector<double> getYvector() const { return Ynew; }

private: // data

  /// vector of the new y-coordinates
  std::vector<double> Ynew;
};

//--------------------------------------------------------------------------//

#endif // MovingAverageFilter_hh

