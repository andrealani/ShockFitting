// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyMeshGenerator_hh
#define ShockFitting_DummyMeshGenerator_hh

//--------------------------------------------------------------------------//

#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyMshGenerator, whose task is read informations
/// needed to generate a mesh from input files and store them in arrays 
/// and vectors.

class DummyMeshGenerator : public MeshGenerator {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DummyMeshGenerator(const std::string& objectName);

  /// Destructor
  virtual ~DummyMeshGenerator();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// Read input files
  virtual void generate();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_MeshGenerator_hh
