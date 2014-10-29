// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_VibrEnergy_hh
#define ShockFitting_VibrEnergy_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include <vector>
#include <string>

//--------------------------------------------------------------------------//

namespace ShockFitting {

/// This class defines a VibrEnergy, whose task is to compute the vibrational
/// energy value.

class VibrEnergy {
public:

  /// Constructor
  VibrEnergy();

  /// Destructor
  ~VibrEnergy();

  /// compute vibrational energy value
  void callVibrEnergy(double, std::vector<double>);

  /// return vibrational energy value
  double getEv() const {return ev;}

private: // helper functions

  /// assign variables used in VibrEnergy to PhysicsData pattern
  void setPhysicsData();

private: //data

  /// vibrational energy
  double ev;

  /// number of species
  unsigned* nsp;

  /// characteristic vibrational temperature (K)
  std::vector<double>* thev;

  /// molecular weight of the species (kg/mol)
  std::vector<double>* mm;

  /// molecules types
  std::vector<std::string>* typemol;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_VibrEnergy_hh
