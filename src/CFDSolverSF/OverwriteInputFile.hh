// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium

// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_OverwriteInputeFile_hh
#define ShockFitting_OverwriteInputFile_hh

//--------------------------------------------------------------------------//

#include<fstream>
#include<sstream>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a OverwriteInputFile, whose task is to change some
/// values inside the CFD solver input file. 
/// The name of values can be specified inside the input.case

class OverwriteInputFile {
public:

  /// Constructor
  /// @param string specifies the name of the in/out files 
  OverwriteInputFile(std::string, std::string);

  /// Destructor
  ~OverwriteInputFile();

  /// overwrite the input file with the new iter value
  /// @param string specifes the name of the variable inside the input case
  /// @param string the replacing value
  void overwriteValue(std::string, std::string);

  /// overwrite two values
  void overwriteValue(std::string, std::string, std::string, std::string);

private: // data

  /// name of the input file
  std::string m_inputFile;

  /// name of the output file
  std::string m_outputFile;

  /// fstream variable reading from input file
  std::ifstream fileIn;

  /// fstream variable writing on output file
  std::ofstream fileOut;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_OverwriteInputeFile_hh 
