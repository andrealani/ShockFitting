// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyCopyMaker_hh
#define ShockFitting_DummyCopyMaker_hh

//--------------------------------------------------------------------------//

#include "Framework/CopyMaker.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyCopyMaker, whose task is to copy mesh
/// arrays

class DummyCopyMaker : public CopyMaker {
public:

  /// Constructor
  /// @param objectName the concrete class name
  DummyCopyMaker (const std::string& objectName);

  /// Destructor
  virtual ~DummyCopyMaker();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// Command arrays copy
  virtual void copy();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_DummyCopyMaker_hh
