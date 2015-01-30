// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium

// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFDSolverSF/OverwriteInputFile.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

OverwriteInputFile::OverwriteInputFile(std::string InputFileName,
                                       std::string OutputFileName)
{
  m_inputFile = InputFileName;
  m_outputFile = OutputFileName;
}

//--------------------------------------------------------------------------//

OverwriteInputFile::~OverwriteInputFile()
{
}

//--------------------------------------------------------------------------//

void OverwriteInputFile::overwriteValue(string VarName, string VarValue)
{
  string dumstring;

  fileIn.open(m_inputFile.c_str());
  fileOut.open(m_outputFile.c_str());

  while(!fileIn.eof()) {
   // read until the VarName is not identified 
   getline(fileIn,dumstring);
   fileOut << dumstring << endl;
   if(dumstring==VarName) {
    // read the old value 
    getline(fileIn,dumstring);
    // write the new value
    fileOut << VarValue << endl;
   }
  }

  fileIn.close();
  fileOut.close();
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
