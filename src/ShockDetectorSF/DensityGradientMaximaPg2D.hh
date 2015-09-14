// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DensityGradientMaximaPg2D_hh
#define ShockFitting_DensityGradientMaximaPg2D_hh

//--------------------------------------------------------------------------//

#include <iostream>
#include <cmath>

//--------------------------------------------------------------------------//

#include "Framework/ShockDetector.hh"
#include "Framework/FileLogManip.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DensityGradientMaximaPg2D, whose task is to find the
/// location of shock discontinuities within a CFD solution as proposed by
/// Pagendarm and Seitz
/// (ref. "An algorithm for detection and visualization of discontinuities
///       in scientific data field applied to flow data with shock waves"
///       H.G. Pagendarm, B. Seitz)

/// @author Valentina De Amicis

class DensityGradientMaximaPg2D : public ShockDetector {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DensityGradientMaximaPg2D (const std::string& objectName);

  /// Destructor
  virtual ~DensityGradientMaximaPg2D();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command finding discontinuities 
  virtual void detect();

  /// command finding discontinuities with input/output arguments
  /// @param primitiveVar  primitive variables in the mesh points  
  virtual void detect(std::vector<double>& primitiveVar);

private: // helper functions

  /// return class name
  std::string getClassName() {return std::string("DensityGradientMaximaPg2D");}

  /// set pressureID
  void setPressureID() {pID=0;}

  /// set the velocityID
  void setVelocityID() {uID=1; vID=2;}

  /// set temperatureID
  void setTemperatureID() {tID=3;}

  /// resize vectors and array
  void setSize();

  /// assign starting pointers for array
  void setAddress();

  /// assign variables used inside DensityGradientMaximaPg2D to MeshData
  void setMeshData();

  /// assign variables used inside DensityGradientMaximaPg2D to PhysicsData
  void setPhysicsData();

  /// de-allocate dynamic array
  void freeArray();

private: // data

  /// value of the filtering factor specified inside the configuration file
  double m_filteringFactor;

  /// x-velocity ID
  unsigned uID;

  /// y-velocity ID
  unsigned vID;

  /// pressure ID
  unsigned pID;

  /// temerature ID
  unsigned tID;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// mesh points ccorinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// working array storing mesh points coordinates
  Array2D<double>* XY;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// vector storing center of each median dual cell
  std::vector<unsigned>* medianDualCellNode;

  /// vector storing pointers to the first surrounding element-node
  std::vector<unsigned>* medianDualCellPtr;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_DensityGradientMaximaPg2D_hh
