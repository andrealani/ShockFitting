// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_RdDpsEq_hh
#define ShockFitting_RdDpsEq_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a RdDpsEq, whose task is to re-distribute shock points 
/// in each shock with uniform spacing

class RdDpsEq : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  RdDpsEq(const std::string& objectName);

  /// Destructor
  virtual ~RdDpsEq();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// re-distribute shock points
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "RdDpsEq";}

  ///assign PhysicsData values to RdDpsEq
  void setPhysicsData();

  ///assign MeshData values to RdDpsEq
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

  /// compute length of shock edge
  void computeShEdgeLength(unsigned); 

  /// compute number of shock points of redistribution
  /// and check the new number of shock points
  void computeNbDistrShPoints(unsigned);

  /// compute distribution step
  void computeDistrStep(unsigned);

  /// interpolation of new shock points
  /// set the first and the last point
  void newPointsInterp(unsigned);

  /// computation of internal points
  void computeInterPoints(unsigned);

  /// rewrite the state arrays and indices
  void rewriteValues(unsigned);

private: // data

  /// dummy variable used as index
  unsigned ISH;

  /// shock edge length
  double Sh_Edge_length;

  /// number of shock points of redistribution
  unsigned nShockPoints_new;

  /// array stores new shock points 
  Array2D <double> XYSh_New;

  /// array stores new zroe upstream values
  Array2D <double> ZRoeShu_New;

  /// array stores new zroe downstream values
  Array2D <double> ZRoeShd_New;

  /// shock length edge 
  std::vector <double> Sh_ABSC;

  /// new shock length edge
  std::vector <double> Sh_ABSC_New;

  /// space dimension
  unsigned* ndim;

  /// nof degree of freedom
  unsigned* ndof;

  /// max nof degree of freedom
  unsigned* ndofmax;

  /// max number of shocks
  /// (needed to set address)
  unsigned* nshmax;

  /// max nof points for each shock
  unsigned* npshmax;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// length of the shock edges in input.case
  double* dxcell;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// upstream state
  Array3D <double>* ZRoeShu;

  /// downstream state
  Array3D <double>* ZRoeShd;

  /// mesh points status
  std::vector <double>* zroeVect;

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_RdDpsEq_hh

