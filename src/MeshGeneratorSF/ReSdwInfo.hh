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
#include "Framework/FileLogManip.hh"
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

  /// read shock info file
  virtual void generate();

  /// Read a given file
  virtual void generate(std::string);

private: // helper functions

  /// get class name
  std::string getClassName() const {return "ReSdwInfo";}

  /// get the input file name
  std::string getInputFiles() const;

  /// read shock infos
  void readShockInfo();

  /// set array characterizing special points
  void setSHinSPPs(unsigned, unsigned);

  /// assign values read and used by ReSdwInfo to PhysicsData
  void setPhysicsData();

  /// assign MeshData values to ReSdwInfo
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

private: //data

  /// number of degree of freedom
  unsigned* ndof;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// type of shock
  std::vector <std::string>* typeSh;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// code characterizing shock points
  Array2D <int>* NodCodSh;

  /// upstream status
  Array3D <double>* ZRoeShu;

  /// downstream status
  Array3D <double>* ZRoeShd;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// old upstream status
  Array3D <double>* ZRoeShuOld;

  /// old downstream status
  Array3D <double>* ZRoeShdOld;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// code characterizing mesh points 
  std::vector <int>* nodcod;

  /// mesh points status
  std::vector <double>* zroe;

  /// reading file variable
  std::ifstream file;

  ///store information on the log file 
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_ReSdwInfo_hh
