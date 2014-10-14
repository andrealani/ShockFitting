// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MeshGenerator_hh
#define ShockFitting_MeshGenerator_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/BaseShockFitting.hh"
#include "Framework/Field.hh"
#include "SConfig/SharedPtr.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

///This class defines a MeshGenerator, whose task is read informations
/// needed to generate a mesh from input files and store them in arrays 
/// and vectors.

class MeshGenerator : public BaseShockFitting {
public:

  ///typedef needed by the self-registration mechanism
  typedef SConfig::Provider<MeshGenerator> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  MeshGenerator(const std::string& objectName);

  /// Destructor
  virtual ~MeshGenerator();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up this object after its last use
  virtual void unsetup() = 0;

  /// Configure the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfgArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// Read one or more input files and stores values in arrays
  virtual void generate() = 0;

  /// Set the mesh field
  virtual void setMeshField (Field *const field);

  /// Get the mesh field
  virtual void getMeshField (Field* field);

  /// Gets the Class name
  static std::string getClassName() {return "MeshGenerator";}

protected: //data

  /// name of the input file
  std::vector<std::string> m_inputFile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_MeshGenerator_hh

