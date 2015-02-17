// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeStateDps_hh
#define ShockFitting_ComputeStateDps_hh

//----------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a ComputeState, whose task is to enforce
/// Rankine-Hugoniot relations to update the solution within the
/// grid-points located on the discontinuities
/// ComputeStatePG     : if NDOF=4 && MODEL="PG"
/// ComputeState4Ar    : if MIXTURE=ar4 && MODEL="Cneq"
/// ComputeState4TCneq : if MODEL=TCneq

class ComputeStateDps : public StateUpdater {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<ComputeStateDps> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  ComputeStateDps(const std::string& objectName);

  /// Destructor
  virtual ~ComputeStateDps();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// update state
  virtual void update() = 0;

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize discontinuity speed array
  void setDiscSpeedSize();

  /// assign values used in ComputeState to MeshData pattern
  void setMeshData();

  /// assign values used in ComputeState to PhysicsData pattern
  void setPhysicsData();

  /// de-allocate the dynamic arrays
  void freeArray();

protected: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of species
  unsigned* nsp;

  /// global indeces
  unsigned* IE;
  unsigned* IEV;
  unsigned* IX;
  unsigned* IY;

  /// heat specific ratio
  double* gref;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;
  
  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* nShockEdges;

  /// type of shock
  std::vector <std::string>* typeSh;

  /// mesh points status
  std::vector<double>* zroeVect;

  /// formation enthalpy at 0K of the species (J/kg)
  std::vector <double>* hf;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// shock/discontinuity speed
  Array3D <double>* WSh;

  /// upstream status
  Array3D <double>* ZroeShu;

  /// downstream status
  Array3D <double>* ZroeShd;

  /// old upstream status
  Array3D <double>* ZroeShuOld;

  /// old downstream status
  Array3D <double>* ZroeShdOld;

  /// total number of shock points
  unsigned TotnbShockPoints;

  /// dummy variables storing normal vector values
  double dx; double dy;

  /// discontinuity speed
  double WS;

  /// working array storing riemann invariant
  Array2D <double>R2;

  /// helping variables
  double help;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

}

#endif // ShockFitting_ComputeStateDps_hh
