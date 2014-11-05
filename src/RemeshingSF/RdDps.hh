// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_RdDps_hh
#define ShockFitting_RdDps_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a RdDps, whose task is to redistribute the shock 
/// points
/// (!) RdDpsEq and RdDps are different classes

class RdDps : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  RdDps(const std::string& objectName);

  /// Destructor
  virtual ~RdDps();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// redistribute shock points
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "RdDps";}

  /// assign variables used in ReDps to MeshData pattern
  void setMeshData();

  /// assign variables used in ReDps to PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointers for arrays 2D and 3D
  void setAddress();

private: // data

  /// space dimension
  unsigned* ndim;

  /// number of degrees of freedom
  unsigned* ndof;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of shock points for each shock
  unsigned* npshmax;

  /// number of shocks
  unsigned* nShocks;

  /// length of the shock edges
  double* dxcell;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges
  std::vector<unsigned>* nShockEdges;

  /// mesh points state
  std::vector<double>* zroeVect;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// shock points upstream state
  Array3D <double>* ZroeShu;

  /// shock points downstream state
  Array3D <double>* ZroeShd;

  /// vector storing shock edges length
  std::vector<double> ShEdgeLgth;

  /// store log file infos
  FileLogManip logfile;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_RdDps_hh
