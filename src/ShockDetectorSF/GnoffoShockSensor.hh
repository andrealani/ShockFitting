// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_GnoffoShockSensor_hh
#define ShockFitting_GnoffoShockSensor_hh

//--------------------------------------------------------------------------//

#include<cmath>

#include "Framework/ShockDetector.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a GnoffoShockSensor, whose task is to find the
/// location of shock discontinuities within a CFD solution
/// as proposed by Gnoffo
/// (ref. "Updates to multidimensional flux reconstruction for
///        Hypersonic simulations on tetrahedral grids", P. Gnoffo, AIAA)

/// @author Valentina De Amicis

class GnoffoShockSensor : public ShockDetector {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  GnoffoShockSensor (const std::string& objectName);

  /// Destructor
  virtual ~GnoffoShockSensor();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command finding discontinuities 
  virtual void detect();

  /// command finding discontinuities with input/output arguments
  /// @param primitiveVar   vector storing primitive variables 
  virtual void detect(std::vector<double>& primitiveVar);

private: // helper functions

  /// return the name of the class
  std::string getClassName() { return "GnoffoShockSensor"; }

  /// minimum value of the pressure ratio 
  unsigned m_minPressRatio;

  /// maximum value of the pressure ratio
  unsigned m_maxPressRatio;

  /// assign varibales used in the sensor to MeshData
  void setMeshData();

  /// assign varibales used in the sensor to PhysicsData
  void setPhysicsData();

  /// assign starting pointers to the Array2D
  void setAddress();

  /// resize vector and arrays
  void setSize();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// specify the fitting strategy
  std::string m_fittingTechnique;

  /// specify the order of the fitting polynomial
  unsigned m_polynomialOrder;

  /// number of vertices for each mesh element
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// IDvariable used for the gnoffo sensor
  unsigned varID;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// specify the type of discontinuity
  std::vector<std::string>* typeSh;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// mesh points ccorinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// working array storing mesh points coordinates
  Array2D<double>* XY;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// shock points coordinates
  Array3D <double>* XYSh;  

  /// array storing the minimum value of the pressure
  /// inside each element and the corresponding node-ID
  Array2D <double> map_varMax;

  /// array storing the maximum value of the pressure
  /// inside each element and the corresponding node-ID
  Array2D <double> map_varMin;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_GnoffoShockSensor_hh
