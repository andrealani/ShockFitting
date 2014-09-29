// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReadPointCloud_hh
#define ShockFitting_ReadPointCloud_hh

//--------------------------------------------------------------------------//

#include "Framework/Connectivity.hh"
#include "Framework/FileProcessing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ReadPointCloud which reads a point cloud from file 
/// and store its data.
/// 
/// @author Andrea Lani
  
class ReadPointCloud : public FileProcessing {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  ReadPointCloud(const std::string& objectName);
  
  /// Destructor
  virtual ~ReadPointCloud();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the file processing 
  virtual void process();
    
  /// Get the solution field
  virtual void getSolutionField(Field* field)
  {
    Connectivity conn; // default connectivity (empty)
    field->reset(m_nbPoints, m_nbVars, conn, &m_array[0]);
  }
    
private: // helper functions
  
  /// get the input file name
  std::string getInputFile() const 
  {
    assert(m_inputFiles.size() == 1); 
    return std::string(m_inputFiles[0] + "." + SConfig::to_str(m_iter));
  }
  
private: // data
  
  /// number of equations
  unsigned m_nbVars;
  
  /// total number of nodes
  unsigned m_nbPoints;

  /// current iteration number
  unsigned m_iter;
  
  /// solution array
  std::vector<double> m_array;
  
  /// number of dummy solution vars to be written to file
  unsigned m_nbDummySolVars;
  
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_ReadPointCloud_hh
