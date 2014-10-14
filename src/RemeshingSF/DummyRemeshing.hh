// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyRemeshing_hh
#define ShockFitting_DummyRemeshing_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyRemeshing, whose task is to alter a mesh
/// geometry

class DummyRemeshing : public Remeshing {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DummyRemeshing (const std::string& objectName);

  /// Destructor
  virtual ~DummyRemeshing();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// command geometry remesh
  virtual void remesh();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_DummyRemeshing_hh
