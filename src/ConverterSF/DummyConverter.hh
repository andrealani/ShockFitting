// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyConverter_hh
#define ShockFitting_DummyConverter_hh

//--------------------------------------------------------------------------//

#include "Framework/Converter.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyConverter, whose task is to convert file format
/// into another format.

class DummyConverter : public Converter {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  DummyConverter(const std::string& objectName);

  /// Destructor
  virtual ~DummyConverter();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// Convert a file format
  virtual void convert();

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_DummyConverter_hh
