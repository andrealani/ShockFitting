// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoNorm_hh
#define Shockfitting_CoNorm_hh

//----------------------------------------------------------------------------//

#include <cmath>
#include <vector>
#include "Framework/Remeshing.hh"
#include "MathTools/Array3D.hh"

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CoNorm, whose task is to compute the normal unit
/// vectors to the shocks and discontinuities in shock/discontinuity points
/// CoNorm4B    : if NDOF=1 && MODEL="B"
/// CoNormPG    : if NDOF=4 && MODEL="PG"
/// CoNorm4Ar   : if MIXTURE=ar4 && MODEL="Cneq"
/// CoNormTCneq : if MODEL=TCneq

class CoNorm : public Remeshing {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<CoNorm> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  CoNorm(const std::string& objectName);

  /// Destructor
  virtual ~CoNorm();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// compute normal vectors
  virtual void remesh() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "CoNorm";}

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

  /// assign values used in CoNorm to MeshData pattern
  void setMeshData();

  /// assign values used in CoNorm to PhysicsData pattern
  void setPhysicsData();

protected: // data

  /// dummy variables for the shock speed
  double ush,vsh;
  double ui, vi;

  /// dummy variables for the shock points coordinates
  double xi,yi,xj,yj,xj2,yj2;

  /// dummy variables for the normal vectors computation
  double taux, tauy, tau;
  double tauxip1, tauyip1, tauxip2, tauyip2;
  double tauxim1, tauyim1, tauxim2, tauyim2;
  double lp1, lp2, lm1, lm2, lp12, lm12, lm22, lp22;
  double dum, nx1, nx2, ny1, ny2, nx4, ny4;

  /// dummy variables for the status recovering
  double uj, vj, roj, help, aj, pj;

  /// dummy variables for the shock dependencies
  int depip1, depim1;

  /// dummy variables for shock indeces
  unsigned I;
  std::vector<unsigned> ISH;
  std::vector<unsigned> IP;

  /// heat specific ratio
  double* gam;

  /// gm1 = gam-1
  double* gm1;

  /// space dimension
  unsigned* ndim;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of species
  unsigned* nsp;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// number of mesh points
  unsigned* npoin;

  /// global indeces
  unsigned* ie;
  unsigned* iev;
  unsigned* ix;
  unsigned* iy;

  /// heat specific ratio
  double* gref;

  /// gas model
  std::vector<std::string>* model;

  /// gas mixture
  std::vector<std::string>* mixture;

  /// number of shocks
  unsigned* r_nShocks;

  /// number of special points
  unsigned* r_nSpecPoints;

  /// number of shock points for each shock
  std::vector<unsigned>* r_nShockPoints;

  /// type of shock
  std::vector<std::string>* r_typeSh;

  /// type of special points
  std::vector<std::string>* r_typeSpecPoints;

  /// mesh points status
  std::vector <double>* zroe;

  /// formation enthalpy at 0K of the species (J/kg)
  std::vector <double>* hf;

  /// constant gas for each species
  std::vector <double>* Rs;

  /// species heat specific ratios
  std::vector <double>* gams;

  /// shock points coordinates
  Array3D <double>* r_XYSh;

  /// downstream status
  Array3D <double>* r_ZRoeShd;

  /// shock points normal vectors
  Array3D <double>* r_vShNor;

  /// array characterizing special points
  Array3D <unsigned>* r_SHinSPPs;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoNorm_hh
