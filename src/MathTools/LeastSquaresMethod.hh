// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef LeastSquaresMethod_hh
#define LeastSquaresMethod_hh

//--------------------------------------------------------------------------//

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "MathTools/Solg.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

/// This class defines a LeastSquaresMethod, whose task is to determine
/// the best fit line to given data
/// up-to-date:
/// "fitEllipse" method is defined to interpolate a point
/// distribution that fit an ellipse shape
/// "fitData" method is defined to fit a points distribution using 
/// polynomials
/// "fitSplittingCurves" is defined to fit a points distribution. In order
/// to accomplish that, the curve is splitted in smaller parts and 
/// different polynomial order can be applied
 
/// @author Valentina De Amicis

class LeastSquaresMethod {
public:

  /// Constructor
  LeastSquaresMethod();

  /// Destructor
  ~LeastSquaresMethod();

  /// find coefficients which best fit the given data
  /// @param X,Y      coordinates of the points distribution
  /// @param polGrad  order of the polynomial
  void fitData(std::vector<double> X,
               std::vector<double> Y,
               unsigned polGrad);

  /// compute coefficients which best fit an ellipse
  /// as referenced in
  /// "Fitting Ellipses and General Second-Order Curves
  /// G. J. Agin, pag.8-10,1981
  /// @param X,Y      coordinates of the points distribution
  void fitEllipse(std::vector<double> X,
                  std::vector<double> Y);

  /// interpolate the given data by splitting the point distribution
  /// in several curves and fitting them with the least squares method 
  /// @param X,Y            coordinates of the points distribution
  /// @param nbSegments     vector of size equal to the physical dimension
  ///        nbSegments(0)  stores the number of x-segments in which the
  ///                       curve will be divided into
  ///        nbSegments(1)  stores the number of y-segments in which the
  ///                       curve will be divided into
  /// @param SegPolynomialOrders  vector of size equal to nbSegments(0)*nbSegments(1)
  ///                             stores the polynomial order that must be used to
  ///                             interpolate each segment of the splitted curve
  /// @param SmoothingOption  specify if the Y-coordinates are smoothed after the
  ///                          least squares method
  void fitSplittingCurves(std::vector<double> X,
                          std::vector<double> Y,
                          std::vector<unsigned> nbSegments,
                          std::vector<unsigned> SegPolynomialOrders,
                          bool SmoothingOptions);
 
  /// return the new components of x
  std::vector<double> getXvector() const { return Xnew; }

  /// return the new components of y
  std::vector<double> getYvector() const { return Ynew; }

  /// return number of shock points
  unsigned getNewPoints() const { return m_nbPoints; }

private: // helper functions

  /// find coefficients which best fit the given data
  /// @param X,Y      coordinates of the points distribution
  /// @param polGrad  order of the polynomial
  void m_fitData(std::vector<double>& X,
                 std::vector<double>& Y,
                 unsigned polGrad);

private: // data

  /// vector of the new x-coordinates
  std::vector<double> Xnew;
 
  /// vector of the new y-coordinates
  std::vector<double> Ynew;

  /// number of shock points
  unsigned m_nbPoints;

  /// vector of the knows variables
  std::vector<double> b;

  /// matrix of the coefficients
  Array2D<double> a;

  /// vector storing coefficients that best fit the data
  std::vector<double> bestCoefficients;
};

//--------------------------------------------------------------------------//
//

#endif // LeastSquaresMethod_hh

