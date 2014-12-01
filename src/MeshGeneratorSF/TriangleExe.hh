// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_TriangleExe_hh
#define ShockFitting_TriangleExe_hh

//--------------------------------------------------------------------------//

#include <vector>
#include <string>
#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a TriangleExe, whose task is to call Triangle
/// Mesh Generator.

class TriangleExe : public MeshGenerator {
public:

  /// Constructor
  TriangleExe(const std::string& objectName);

  /// Destructor
  virtual ~TriangleExe();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// generate the new mesh
  virtual void generate();

private:

  /// dummy string
  std::string command;

  /// name of current file
  std::stringstream* fname;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting_TriangleExe_hh

//--------------------------------------------------------------------------//

#endif // ShockFitting_TriangleExe_hh
