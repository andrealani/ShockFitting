// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Param2Prim_hh
#define ShockFitting_Param2Prim_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include <vector>
#include <string>

#include "MathTools/Array2D.hh"
#include "Framework/VariableTransformer.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Param2Prim, whose task is to operate a variable
/// transformation from Roe parameter vector variables to primitive variables.
/// According to gas model (Pg, Cneq or TCneq) and to CFmesh output 
/// (dimensional or adimensional) several objects are defined:
/// .) Param2PrimPgAdimensional
/// .) Param2PrimPgDimensional
/// .) Param2PrimTCneqAdimensional
/// .) Param2PrimTCneqDimesnional
/// .) Param2PrimCneqAdimensional
/// .) Param2Prim4CneqDimensional

class Param2Prim : public VariableTransformer {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<Param2Prim> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  Param2Prim(const std::string& objectName);

  /// Destructor
  virtual ~Param2Prim();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// command variables transformation
  virtual void transform() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "Param2Prim";}

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// assign variables used in Parm2Prim to MeshData pattern
  void setMeshData();

  /// assign variables used in Parm2Prim to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

protected: // data

  /// space dimension
  unsigned* ndim;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of mesh points
  unsigned* npoin;

  /// number of species
  unsigned* nsp;

  /// gas constant (value read by ReferenceInfo)
  double* Rgas;

  /// heat specific ratio (value read by ReferenceInfo)
  double* gam;

  /// freestream temperature
  double* Tref;

  /// freestream pressure  
  double* pref;

  /// reference density 
  double* rhoref;

  /// freestream velocity
  double* uref;

  /// reference length
  double* Lref;

  /// global indeces
  unsigned* ix;
  unsigned* iy;
  unsigned* ie;
  unsigned* iev;

  /// number of molecules
  unsigned* nmol;

  /// molecular weight of the species (kg/mol)
  std::vector <double>* mm;

  /// formation enthalpy at 0K of the species (J/kg)
  std::vector <double>* hf;

  /// characteristic vibrational temperature (K)
  std::vector <double>* thev;

  /// specific heat ratio for each species
  std::vector <double>* gams;

  /// molecules types
  std::vector <std::string>* typemol;

  /// mesh points status
  std::vector<double>* zroe;

  /// mesh points coordinates
  std::vector<double>* coor;

  /// mesh points status (in array storing)
  Array2D <double>* v_Zroe;

  /// mesh points coordinates (in array storing)
  Array2D <double>* v_XY;

  /// dummy variables used to make the transformation
  double rho; double kinetic; double h; double help;
  double pres; double p; double u; double v;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Param2Prim_hh
