// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_RedistributeShockPoints_hh
#define ShockFitting_RedistributeShockPoints_hh

//--------------------------------------------------------------------------//

#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a RedistributeShockPoints, whose task is to 
/// re-distribute shock points in each shock with uniform spacing

class RedistributeShockPoints {
public:

  /// Constructor
  RedistributeShockPoints();

  /// Destructor
  ~RedistributeShockPoints();

  /// re-distribute shock points
  void distribute();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "RedistributeShockPoints";}

  ///assign PhysicsData values to RedistributeShockPoints
  void setPhysicsData();

  ///assign MeshData values to RedistributeShockPoints
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

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// dummy variable used as index
  unsigned ISH;

  /// shock edge length
  double Sh_Edge_length;

  /// number of shock points of redistribution
  unsigned nShockPoints_new;

  /// array stores new shock points 
  Array2D <double> XYSh_New;

  /// shock length edge 
  std::vector <double> Sh_ABSC;

  /// new shock length edge
  std::vector <double> Sh_ABSC_New;

  /// nof degree of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_RedistributeShockPoints_hh

