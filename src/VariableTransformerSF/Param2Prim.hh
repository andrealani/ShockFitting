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
/// The new values are stored in new arrays to not overwrite the old mesh status.
/// These new values are pushed back at the end of the arrays of the old mesh
/// in the fortran version these new arrays are referred to index "1"
/// (ex: ZROE(1))
///
/// According to the gas model (Pg, Cneq or TCneq) and to CFmesh output 
/// (dimensional or adimensional) several objects are defined:
/// .) Param2PrimPgAdimensional
/// .) Param2PrimPgDimensional
/// .) Param2PrimTCneqAdimensional (not implemented yet)
/// .) Param2PrimTCneqDimesnional
/// .) Param2PrimCneqAdimensional (not implemented yet)
/// .) Param2PrimCneqDimensional (not implemented yet)

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

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of shock points for each shock
  unsigned* npshmax;

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
  unsigned* IX;
  unsigned* IY;
  unsigned* IE;
  unsigned* IEV;

  /// number of molecules
  unsigned* nmol;

  /// number of mesh points
  std::vector<unsigned>* npoin;

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
  std::vector<double>* zroeVect;

  /// mesh points coordinates
  std::vector<double>* coorVect;

  /// mesh points status (in array storing)
  Array2D <double>* zroe;

  /// mesh points coordinates (in array storing)
  Array2D <double>* XY;

  /// dummy variable to store new arrays size
  unsigned totsize;

  /// dummy variable to store array starting pointer
  unsigned start;

  /// dummy variables used to make the transformation
  double rho, kinetic, h, help;
  double pres, p, u, v;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Param2Prim_hh
