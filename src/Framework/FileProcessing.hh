// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FileProcessing_hh
#define ShockFitting_FileProcessing_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/BaseShockFitting.hh"
#include "Framework/Field.hh"
#include "SConfig/SharedPtr.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {
  
//--------------------------------------------------------------------------//

/// This class defines a FileProcessing, whose task is to transform all sort 
/// of data file manipulations, including reordering, conversion from ASCII to
/// binary (and viceversa), translation from one datastructure to another, 
/// modifications of headers, file joining or splitting, etc. 
/// 
/// @author Andrea Lani

class FileProcessing : public BaseShockFitting {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<FileProcessing> PROVIDER;
  
  /// Constructor 
  /// @param objectName the concrete class name
  FileProcessing(const std::string& objectName);
  
  /// Destructor
  virtual ~FileProcessing();
   
  /// Set up this object before its first use
  virtual void setup() = 0;
  
  /// Unset up this object after its last use
  virtual void unsetup() = 0;
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
  /// Process one or more input files and save the result into one or more output files
  virtual void process() = 0;
  
  /// Get the solution field
  virtual void getSolutionField(Field* field);
  
  /// Set the solution field
  virtual void setSolutionField(Field *const field);
  
  /// Get the mesh field (coordinates)
  virtual void getMeshField(Field* field);
  
  /// Set the mesh field
  virtual void setMeshField(Field *const field);
  
  /// Gets the Class name
  static std::string getClassName() {return "FileProcessing";}
  
protected: // functions
  
  /// Rescale a given field
  void rescaleField(Field *const field);
  
  /// Get starting iteration
  unsigned getStartIter() const {return m_startIter;}

  /// Get starting iteration
  unsigned getIORate() const {return m_ioRate;}
  
protected: // data
  
  /// names of the input files
  std::vector<std::string> m_inputFiles;
  
  /// names of the output files
  std::vector<std::string> m_outputFiles;

  /// names of the I/O variables
  std::vector<std::string> m_varNames;

  /// scaling factors for each I/O variable
  std::vector<double> m_varScalingFactors;
  
  /// flag telling if to I/O dummy data
  bool m_dummyData;

  /// starting iteration
  unsigned m_startIter;

  /// rate telling how often a file is processed
  unsigned m_ioRate; 
  
};

//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_FileProcessing_hh
