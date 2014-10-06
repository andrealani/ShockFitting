// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ChemConst_hh
#define ShockFitting_ChemConst_hh

//--------------------------------------------------------------------------//

#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ChemConst, whose tasks is to read chemical parameters
/// and assign them to PhysicsData pattern.

class ChemConst: public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ChemConst(const std::string& objectName);

  /// Destructor
  virtual ~ChemConst();

  /// Set up object before its first use
  virtual void setup();

  /// Unset up object before its last use
  virtual void unsetup();

  /// Run constant values reading
  virtual void generate();

private: //helper functions

  /// assign values read by ChemParam to PhysicsData
  void setPhysicsData();

private: //data

  /// gas constant [J*mole-1*Ke-1]
  double m_R;

  /// avogadro constant
  double m_Na;

  /// Boltzmann constant [eV*Ke-1]
  double m_K;

  /// Boltzmann constant [J*Ke-1]
  double m_KSI;

  ///
  double* r_R;

  /// 
  double* r_Na;

  ///
  double* r_K;

  ///
  double* r_KSI;
};

} // namespace Shockfitting

#endif //ShockFitting_ChemConst_hh
