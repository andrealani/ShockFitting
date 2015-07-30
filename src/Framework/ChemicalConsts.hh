// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ChemicalConsts_hh
#define ShockFitting_ChemicalConsts_hh

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// Provides an a set of static functions for chemical constants.

class ChemicalConsts {
public:

  /// Universal constant of gas [J/mol/K]
  static double Rgp() {return 8.314472;}

  /// Avogadro Constant [1/mol]
  static double Na() {return 6.0221367e23;}

  /// Boltzmann constant [ev/K]
  static double KBol() {return 8.617386e-5;}
  
  /// Boltzmann constant SI [J/K]
  static double KBolSI() {return 1.380651488e-23;}
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // Shockfitting_ChemicalConsts
