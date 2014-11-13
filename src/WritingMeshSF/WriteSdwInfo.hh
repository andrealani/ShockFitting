// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteSdwInfo_hh
#define ShockFitting_WriteSdwInfo_hh

//--------------------------------------------------------------------------//

#include <fstream>
#include <vector>
#include "Framework/FileLogManip.hh"
#include "Framework/WritingMesh.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a WriteSdwInfo, whose task is to write file
/// sh99.dat containg informations about shocks and or discontinuities

class WriteSdwInfo : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteSdwInfo(const std::string& objectName);

  /// Destructor
  virtual ~WriteSdwInfo();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// write the sh99 file
  virtual void write();

private: // helper functions

  /// return class name
  std::string getClassName() const { return std::string("WriteSdwInfo"); }

  /// assign variables used in WriteSdwInfo to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteSdwInfo to PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointers to Array3D
  void setAddress();

  /// write SHinSPPs elements of sh99 file
  void writeSHinSPPs(unsigned, unsigned);

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number oer shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* nShockEdges;

  /// type of shock
  std::vector<std::string>* typeSh;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// mesh points state
  std::vector<double>* zroe;

  /// shock points upstream state
  Array3D <double>* ZroeShu;

  /// shock points downstream state
  Array3D <double>* ZroeShd;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// fstream variable writing on sh99 file
  std::ofstream file;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteSdwInfo.hh
