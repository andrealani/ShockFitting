// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WritePointCloud_hh
#define ShockFitting_WritePointCloud_hh

//--------------------------------------------------------------------------//

#include "Framework/Connectivity.hh"
#include "Framework/FileProcessing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//
  
/// This class defines a WritePointCloud which writes a point cloud to file 
/// from given data
/// 
/// @author Andrea Lani

class WritePointCloud : public FileProcessing {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  WritePointCloud(const std::string& objectName);
  
  /// Destructor
  virtual ~WritePointCloud();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the file processing 
  virtual void process();
    
  /// Get the solution field
  virtual void setSolutionField(Field *const field) {assert(field != NULL); m_field = field;}
  
private: // helper functions
  
  /// get the input file name
  std::string getOutputFile() const 
  {
    assert(m_outputFiles.size() == 1); 
    return std::string(m_outputFiles[0] + "." + SConfig::to_str(m_iter));
  }
  
private: // data
  
  /// number of equations
  Field* m_field;
  
  /// iteration number
  unsigned m_iter;

  /// write file in TECPLOT format
  bool m_writeTecplot;
  
  /// number of points in I,J,K
  std::vector<unsigned> m_nijk;
  
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_WritePointCloud_hh
