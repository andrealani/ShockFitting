// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReadInterpolateWrite_hh
#define ShockFitting_ReadInterpolateWrite_hh

//--------------------------------------------------------------------------//

#include "Framework/ShockFittingObj.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {
  
//--------------------------------------------------------------------------//

/// This class defines a ReadInterpolateWrite, whose task is to read a file, 
/// interpolate a field and write another file with the interpolated field
/// 
/// @author Andrea Lani

class ReadInterpolateWrite : public ShockFittingObj {
public:
    
  /// Constructor 
  /// @param objectName the concrete class name
  ReadInterpolateWrite(const std::string& objectName);
  
  /// Destructor
  virtual ~ReadInterpolateWrite();
  
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the coupling tools
  virtual void process();
  
protected: 

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
private:
  
  /// command object reading the first input file
  SConfig::SharedPtr<FileProcessing> m_readInputFile1; 
  
  /// command object reading the second input file
  SConfig::SharedPtr<FileProcessing> m_readInputFile2; 
  
  /// command object writing the output file
  SConfig::SharedPtr<FileProcessing> m_writeOutputFile; 
  
  /// command object interpolating from input to output
  SConfig::SharedPtr<FieldInterpolator> m_ioInterpolator; 
    
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ReadInterpolateWrite_hh
