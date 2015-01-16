// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Prim2ParamTCneqDimensional_hh
#define ShockFitting_Prim2ParamTCneqDimensional_hh

//--------------------------------------------------------------------------//

#include "VariableTransformerSF/Prim2Param.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Prim2ParamTCneqDimensional, whose task is to 
/// transform primitive variables in dimensional Roe parameter vector
/// variables for a TCneq model.
/// The new values are stored in new arrays to not overwrite the old mesh status.
/// These new values are pushed back at the end of the arrays of the old mesh
/// in the fortran version these new arrays are referred to index "1"
/// (ex: ZROE(1))

class Prim2ParamTCneqDimensional : public Prim2Param {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Prim2ParamTCneqDimensional(const std::string& objectName);

  /// Destructor
  virtual ~Prim2ParamTCneqDimensional();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command variables transformation
  virtual void transform();

  /// command given variables transformation
  virtual void transform(std::vector<double>*, 
                         std::vector<double>*, std::vector<double>*);

private: // data

  /// concentrations of the chemical species
    std::vector<double> alpha;

  /// species densities
  std::vector<double> rhos;

  /// speed components
  std::vector<double> u;

  /// temperatures
  std::vector<double> T;

  /// vibrational energy
  double ev;

  /// dummy index
  unsigned IM;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Prim2ParamTCneqDimensional_hh
