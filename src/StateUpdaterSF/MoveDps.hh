// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MoveDps_hh
#define ShockFitting_MoveDps_hh

//----------------------------------------------------------------------------//

#include "Framework/StateUpdater.hh"
#include "MathTools/Array3D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a MoveDps, whose task is to compute the new 
/// position of the shock front by displacing all discontinuity points
/// using the first order integration formula
/// P(t+dt) = P(t) + w(t+dt)*n*dt
/// MoveDpsPg    : if NDOF=4 && MODEL=Pg
/// MoveDpsAr    : if MODEL = Cneq && MIXTURE = Ar
/// MoveDpsTCneq : if MODEL=TCneq

class MoveDps : public StateUpdater {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<MoveDps> PROVIDER;

  /// Constructor
  MoveDps(const std::string& objectName);

  /// Destructor
  virtual ~MoveDps();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// move shock
  virtual void update() = 0;

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// assign values used in MoveDps to MeshData pattern
  void setMeshData();

  /// assign values used in MoveDps to PhysicsData pattern
  void setPhysicsData();

protected: // data

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
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_MoveDps_hh
