// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoNorm_hh
#define ShockFitting_CoNorm_hh

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

  /// de-allocate dynamic arrays
  void freeArray();

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

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// type of shock
  std::vector<std::string>* typeSh;

  /// type of special points
  std::vector<std::string>* typeSpecPoints;

  /// mesh points status
  std::vector <double>* zroe;

  /// formation enthalpy at 0K of the species (J/kg)
  std::vector <double>* hf;

  /// constant gas for each species
  std::vector <double>* Rs;

  /// species heat specific ratios
  std::vector <double>* gams;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// downstream status
  Array3D <double>* ZRoeShd;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CoNorm_hh
