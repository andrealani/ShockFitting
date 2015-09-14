// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_NormalMachNumberPg2D_hh
#define ShockFitting_NormalMachNumberPg2D_hh

//----------------------------------------------------------------------------//

#include <iostream>
#include <cmath>

//----------------------------------------------------------------------------//

#include "Framework/ShockDetector.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/FileLogManip.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "Framework/VariableTransformer.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a NormalMachNUmber shock detector, whose task is to
/// find the location of shock discontinuities within a CFD solution as
/// proposed by Lovely and Haimes
/// (ref. "Shock Detection from CFD results", D. Lovely, R. Hames, AIAA)
/// It can be used for the two-dimensional flow with one shock

/// @author Valentina De Amicis

class NormalMachNumberPg2D : public ShockDetector {
public:

  /// Constructor
  /// @param objectName the concrete class name
  NormalMachNumberPg2D(const std::string& objectName);

  /// Destructor
  virtual ~NormalMachNumberPg2D();

  /// Setup object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// find the shock location
  virtual void detect();

  /// find the shock location using given vector infos
  /// @param primitiveVar      primitive variables in the mesh points 
  virtual void detect(std::vector<double>& primitiveVar);

private: // helper function

  /// return class name
  std::string getClassName() {return std::string("NormalMachNumberPg2D");}

  /// 
  bool falseResult(double, double);

  /// set pressureID
  void setPressureID() {pID=0;}

  /// set the velocityID
  void setVelocityID() {uID=1; vID=2;}

  /// set temperatureID
  void setTemperatureID() {tID=3;}

  /// resize vectors and array
  void setSize();

  /// assign starting pointers to array
  void setAddress();

  /// assign variables used in NormalMachNumberPg2D to MeshData
  void setMeshData();

  /// assign variables used in NormalMachNumberPg2D to MeshData
  void setPhysicsData();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of vertices in each element
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* nShockEdges;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// mesh points ccorinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes element (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// working array storing mesh points coordinates
  Array2D<double>* XY;

  /// shock points coordinates
  Array3D<double>* XYSh;

  /// x-velocity ID
  unsigned uID;

  /// y-velocity ID 
  unsigned vID;
  
  /// pressure ID
  unsigned pID;

  /// temerature ID
  unsigned tID;

  /// vector storing center of each median dual cell
  std::vector<unsigned>* medianDualCellNode;

  /// vector storing pointers to the first surrounding element-node
  std::vector<unsigned>* medianDualCellPtr;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_NormalMachNumberPg2D_hh
