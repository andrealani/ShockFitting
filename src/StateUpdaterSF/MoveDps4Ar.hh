// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MoveDpsAr4_hh
#define ShockFitting_MoveDpsAr4_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/FileLogManip.hh"
#include "Framework/StateUpdater.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a MoveDps4Ar, whose task is to compute the new 
/// position of the shock for the Cneq model and the ar4 mixture

class MoveDps4Ar : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  MoveDps4Ar(const std::string& objectName);

  /// Destructor
  ~MoveDps4Ar();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// move the shock points
  void update();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("MoveDpsAr");}

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// assign values used in MoveDps to MeshData pattern
  void setMeshData();

  /// assign values used in MoveDps to PhysicsData pattern
  void setPhysicsData();

  /// de-allocate dynamic arrays
 void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// specific heat ratio (for Cneq and TCneq models)
  double* gref;

  /// number of species
  unsigned* nsp;

  /// global indeces
  unsigned* IX;
  unsigned* IY;
  unsigned* IE;
  unsigned* IEV;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges
  std::vector<unsigned>* nShockEdges;

  /// type of discontinuity
  std::vector<std::string>* typeSh;

  /// mesh points state
  std::vector<double>* zroeVect;

  /// formation enthalpy at 0K of the species (J/kg)
  std::vector<double>* hf;

  /// shock points coordinates
  Array3D <double>* XYSh;
  
  /// shock/discontinuity speed
  Array3D <double>* WSh;

  /// shock points state
  Array3D <double>* ZroeSh;

  /// max dt
  double dt;

  /// work variables
  double ro, a, p, help, dum;
  double WShMod;

  /// maximum value of the shock points speed
  double WShMax;

  // store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_MoveDps4Ar_hh
