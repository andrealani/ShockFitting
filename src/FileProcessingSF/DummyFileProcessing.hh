// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DummyFileProcessing_hh
#define ShockFitting_DummyFileProcessing_hh

//--------------------------------------------------------------------------//

#include "Framework/FileProcessing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a DummyFileProcessing, whose task is to transform all sort 
/// of data file manipulations, including reordering, conversion from ASCII to
/// binary (and viceversa), translation from one datastructure to another, 
/// modifications of headers, file joining or splitting, etc. 
/// 
/// @author Andrea Lani

class DummyFileProcessing : public FileProcessing {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  DummyFileProcessing(const std::string& objectName);
  
  /// Destructor
  virtual ~DummyFileProcessing();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Process one or more input files and save the result into one or more output files
  virtual void process();
  
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_DummyFileProcessing_hh
