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
   if(dumstring==VarName) {
    // write the new value
    fileOut << VarValue << endl;
   }
   else { fileOut << dumstring << endl;}
  }

  fileIn.close();
  fileOut.close();
}

//--------------------------------------------------------------------------//

void OverwriteInputFile::overwriteValue(string VarName1, string VarValue1,
                                        string VarName2, string VarValue2)
{
  string dumstring;
  
  fileIn.open(m_inputFile.c_str());
  fileOut.open(m_outputFile.c_str());

  while(!fileIn.eof()) {
   // read until the VarName1 or VarName2 is not identified
   getline(fileIn,dumstring);

   if(dumstring==VarName1) {
    // write the new value
    fileOut << VarValue1;
   }

   else if(dumstring==VarName2) {
    // write the new value
    fileOut << VarValue2;
   }

   else {
    fileOut << dumstring << endl;
   }
  }

  fileIn.close();
  fileOut.close();

}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
