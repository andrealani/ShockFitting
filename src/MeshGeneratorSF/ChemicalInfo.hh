// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ChemicalInfo_hh
#define ShockFitting_ChemicalInfo_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "SConfig/StringManip.hh"
#include "Common/FileLogManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ChemicalInfo, whose task is to take informations about
/// gas model and mixture.
/// If model = Cneq or model = TCneq ChemInfo reads data from mixture.dat
/// file and stores them in vectors 

class ChemicalInfo : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ChemicalInfo(const std::string& objectName);

  /// Destructor
  virtual ~ChemicalInfo();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run the mixture file reading
  virtual void generate();

private: // helper functions

  /// get class name
  std::string getClassName() const {return "ChemicalInfo";}

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

  /// reference speed
  /// used only if MODEL = Cneq && MIXTURE= ar4
  double m_Qref;

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

  /// number of the chemical species
  unsigned* nsp;

  /// name of the species IUPAC
  std::vector <std::string>* name;

  /// molecular weight of the species (kg/mol)
  std::vector <double>* mm;

  /// formation enthalpy at 0K of the species (J/kg
  std::vector <double>* hf;

  /// characteristic vibrational temperature (K)
  std::vector <double>* thev;

  /// specific heat ratio for each species
  std::vector <double>* gams;

  /// type of molecule
  std::vector <std::string>* typemol;

  /// mixture name in the mixture file
  std::string mixtureFile_name;

  /// dummy string
  std::string dummystr;

  /// reading file
  std::ifstream file;

  /// store file log info
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_ChemicalInfo_hh
