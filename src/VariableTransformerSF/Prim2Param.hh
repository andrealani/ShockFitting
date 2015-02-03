// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Prim2Param_hh
#define ShockFitting_Prim2Param_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include <vector>
#include <string>

#include "MathTools/Array2D.hh"
#include "Framework/VariableTransformer.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Prim2Param, whose task is to operate a variable
/// transformation from primitive cariables to Roe parameter vector
/// variables.
/// The new values are stored in new arrays to not overwrite the old mesh status.
/// These new values are pushed back at the end of the arrays of the old mesh
/// in the fortran version these new arrays are referred to index "1"
/// (ex: ZROE(1))
///
/// According to gas model (Pg, Cneq or TCneq) and to Roe Variables 
/// (dimensional or adimensional) several objects are defined:
/// .) Prim2ParamPgAdimesnional
/// .) Prim2ParamPgDimensional
/// .) Prim2ParamTCneqAdimensional (not implemented yet)
/// .) Prim2ParamTCneqDimensional
/// .) Prim2ParamCneqAdimensional (not implemented yet)
/// .) Prim2ParamCneqDimensional (not implemented yet)

class Prim2Param : public VariableTransformer {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<Prim2Param> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  Prim2Param(const std::string& objectName);

  /// Destructor
  virtual ~Prim2Param();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// command variables transformation
  virtual void transform() = 0;

  /// Transform given variables transformation
  virtual void transform(std::vector<double>*, std::vector<double>*,
                         std::vector<double>*) = 0;

  /// Gets the Class name
  static std::string getClassName() {return "Prim2Param";}

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// assign variables used in Prim2Param to MeshData pattern
  void setMeshData();

  /// assign variables used in Prim2Param to PhysicsData pattern
  void setPhysicsData();

  /// assign a few number of variables used in Prim2Param to PhysicsData pattern
  void setPhysicsData(std::string);

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

protected: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of species
  unsigned* nsp;

  /// global indeces
  unsigned* IX;
  unsigned* IY;
  unsigned* IE;
  unsigned* IEV;

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

  /// mesh points status (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// mesh points status (in array storing)
  Array2D <double>* zroe;

  /// mesh points status(in array storing)
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

#endif // ShockFitting_Prim2Param_hh

