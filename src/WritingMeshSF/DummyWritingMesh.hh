// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyWritingMesh_hh
#define ShockFitting_DummyWritingMesh_hh

//--------------------------------------------------------------------------//

#include "Framework/WritingMesh.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyWritingMesh, whose task is to write data
/// computed by the code in a mesh generator format.

class DummyWritingMesh : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  DummyWritingMesh(const std::string& objectName);

  /// Destructor
  virtual ~DummyWritingMesh();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// write a mesh generator file
  virtual void write();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_DummyWritingMesh_hh
