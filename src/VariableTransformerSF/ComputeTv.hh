// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeTv_hh
#define ShockFitting_ComputeTv_hh

//--------------------------------------------------------------------------//

#include "VariableTransformerSF/Param2Prim.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeTv, whose task is to compute the vibrational
/// temperature

class ComputeTv {
public:

  /// Constructor
  ComputeTv();

  /// Destructor
  ~ComputeTv();

  /// compute T values
  void callComputeTv(double, std::vector<double>,
                      std::vector<double>);

  /// return T values
  std::vector<double> getT() const { return T; }

private: // helper functions

  void funct();

  /// assign variables used in ComputeTv to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// vibrational energy
  double ev;
 
  /// concentrations of the chemical species
  std::vector<double> alpha; 

  /// temperature                      
  std::vector<double> T;

  /// number of species
  unsigned* nsp;

  /// characteristic vibrational temperature (K)
  std::vector<double>* thev;
  
  /// molecular weight of the species (kg/mol)
  std::vector<double>* mm;

  /// molecules types
  std::vector<std::string>* typemol;
  
  double F, Fd;
  double tollerance;
  double tolleranceD;
  unsigned KMAX;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeTv_hh
