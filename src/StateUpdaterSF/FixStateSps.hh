// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FixStateSps_hh
#define ShockFitting_FixStateSps_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/StateUpdater.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FixStateSps, whose task is to fix and correct 
/// the  nodal values and shock speed in all special points using
/// the correct S-S interaction relations

class FixStateSps : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  FixStateSps(const std::string& objectName);

  /// Destructor
  virtual ~FixStateSps();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// fix the state in the special points
  virtual void update();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "FixStateSps";}

  /// fix the states for IPX and IPY special points
  void fixIPXandIPYspecPoints(unsigned);

  /// fix the states for OPX and OPY special points
  void fixOPXandOPYspecPoints(unsigned);

  /// fix the states for WPNRX and WPNRY special points
  void fixWPNRXandWPNRYspecPoints(unsigned);

  /// fix the states for TP special point
  void fixTPspecPoints(unsigned);

  /// fix the states for QP special point
  void fixQPspecPoints(unsigned);

  /// fix the states for RRX special point
  void fixRRXspecPoints(unsigned);

  /// fix the states for EP special point
  void fixEPspecPoints(unsigned);

  /// fix the states for C special point
  void fixCspecPoints(unsigned);

  /// set states
  void setState(std::string, unsigned);

  /// set shock slopes
  void setShockSlope(unsigned);

  /// set shock characterizing indeces
  void setShockIndeces(unsigned, unsigned);

  /// assign variables used in FixStateSps to MeshData pattern
  void setMeshData();

  /// assign variables used in FixStateSps to PhysicsData pattern
  void setPhysicsData();

  /// set starting pointers for the arrays
  void setAddress();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
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
  std::vector <std::string>* typeSh;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// mesh points status
  std::vector<double>* zroeVect;

  /// mesh points coordinates
  std::vector<double>* coorVect;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// upstream shock points coordinates
  Array3D <double>* XYShu;

  /// downstream shock points coordinates
  Array3D <double>* XYShd;

  /// upstream status
  Array3D <double>* ZroeShu;

  /// downstream status
  Array3D <double>* ZroeShd;

  /// old upstream status
  Array3D <double>* ZroeShuOld;

  /// old downstream status
  Array3D <double>* ZroeShdOld;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// shock points speed
  Array3D <double>* WSh;

  /// shock indeces
  std::vector<unsigned> ISH;
  std::vector<unsigned> IP;
  unsigned I;

  /// Riemann invariants
  double R23, R14;

  /// work variables storing shock speed
  double WS;
  double WWS;
  double dx, dy;
  double DXR14, DYR14;

  /// work variable setting the start index in the state recovering
  unsigned startIndex;

  /// work variable saving current shock indeces
  unsigned ip, ish;

  /// work variables storing the states
  std::vector<double> xi;
  std::vector<double> x;

  /// store log files infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_FixStateSps_hh
