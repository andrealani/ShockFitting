// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ChemInfo_hh
#define ShockFitting_ChemInfo_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ChemInfo, whose task is to take informations about
/// gas model and mixture.
/// If model = Cneq or model = TCneq ChemInfo reads data from mixture.dat
/// file and stores them in vectors 

class ChemInfo : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ChemInfo(const std::string& objectName);

  /// Destructor
  virtual ~ChemInfo();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run the mixture file reading
  virtual void generate();

private: // helper functions

  /// get the input file name
  std::string getInputFiles() const;

  /// assign values read by ReadTriangle to PhysicsData
  void setPhysicsData();

  /// resize vectors
  void setSize();

  /// set gas model
  void setModel();

  /// set mixture name in the mixture file
  void setMixtureFileName();

  /// set gas mixture
  void setMixture();

  /// set global indeces
  void setGlobalIndex();

  /// set number of degree of freedom
  void setnbDof();

  /// set number of species
  void setnbSpecies();

  /// set names of the species
  void setSpeciesNames();

  /// set species molecular weights
  void setMolWeights();

  /// set formation enthalpy
  void setFormEnthalp();
 
  /// set Vibrational temperatures
  void setVibrTemp();

  /// set specific heat ratios
  void setSpecHeatRatio();

  /// set molecular types
  void setMolTypes();

private: // data

  /// global index
  unsigned m_IE;

  /// global index
  unsigned m_IX;

  /// global index
  unsigned m_IY;

  /// global index
  unsigned m_IEV;

  /// model
  std::string m_model;

  /// mixture
  std::string m_mixture;

  /// space dimension
  unsigned* ndim;

  /// number of degree of freedom
  unsigned* ndof;

  /// global index (assignable to PhysicsData)
  unsigned* ie;

  /// global index (assignable to PhysicsData)
  unsigned* ix;

  /// global index (assignable to PhysicsData)
  unsigned* iy;

  /// global index (assignable to PhysicsData) 
  unsigned* iev;

  ///  model (assignable to PhysicsData)
  std::vector <std::string>* model;

  /// mixture (assignable to PhysicsData)
  std::vector <std::string>* mixture;

  /// number of species
  unsigned* nsp;

  /// name of species
  std::vector <std::string>* name;

  /// species molecular weights [kg*mole-1]
  std::vector <double>* mm;

  /// species formation enthalpy [J*kge-1]
  std::vector <double>* hf;

  /// species vibrational temperatures [K]
  std::vector <double>* thev;

  /// species specif height ratios
  std::vector <double>* gams;

  /// species molecular types
  std::vector <std::string>* typemol;

  /// mixture name in the mixture file
  std::string mixtureFile_name;

  /// reading file
  std::ifstream file;

};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_MixtureInfo_hh
