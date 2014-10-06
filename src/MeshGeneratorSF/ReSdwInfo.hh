// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReSdwInfo_hh
#define ShockFitting_ReSdwInfo_hh

//--------------------------------------------------------------------------//

#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/MeshGenerator.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ReSdwInfo, whose task is read data fromsh00.dat file
/// and store them in vectors and array
///

class ReSdwInfo : public MeshGenerator {
public:

  /// Constructor
  ReSdwInfo(const std::string& objectName);

  /// Destructor
  virtual ~ReSdwInfo();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run file reading
  virtual void generate();

private: // helper functions

  /// get the input file name
  std::string getInputFiles() const;

  /// read shock infos
  void readShockInfo();

  /// assign values read and used by ReSdwInfo to PhysicsData
  void setPhysicsData();

  /// assign MeshData values to ReSdwInfo
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

private: //data

  /// space dimension
  unsigned* ndim;

  /// number of degree of freedom
  unsigned* ndof;

  /// max number of degree of freedom
  unsigned* ndofmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// max number of special points
  unsigned* nspmax;

  /// number of mesh points
  unsigned* npoin;

  /// number of shocks
  unsigned* r_nShocks;

  /// number of special points
  unsigned* r_nSpecPoints;

  /// number of shock points for each shock
  std::vector <unsigned>* r_nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* r_nShockEdges;

  /// type of shock
  std::vector <std::string>* r_typeSh;

  /// type of special points
  std::vector <std::string>* r_typeSpecPoints;

  /// code characterizing shock points
  Array2D <int>* r_NodCodSh;

  /// upstream state
  Array3D <double>* r_ZRoeShu;

  /// downstream state
  Array3D <double>* r_ZRoeShd;

  /// shock points coordinates
  Array3D <double>* r_XYSh;

  /// old upstream state
  Array3D <double>* r_ZRoeShuOld;

  /// old downstream state
  Array3D <double>* r_ZRoeShdOld;

  /// array characterizing special point
  Array3D <unsigned>* r_SHinSPPs;

  /// code characterizing mesh points 
  std::vector <int>* nodcod;

  /// mesh points status
  std::vector <double>* zroe;

  /// reading file variable
  std::ifstream file;

};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_ReSdwInfo_hh
