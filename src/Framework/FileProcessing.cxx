// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FileProcessing.hh"
#include "Framework/Log.hh"
#include "Framework/IOFunctions.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//
  
FileProcessing::FileProcessing(const std::string& objectName) :
  BaseShockFitting(objectName)
{
  m_inputFiles = vector<string>();
  addOption("InputFiles",&m_inputFiles, 
	    "List of the names of input files");  
  
  m_outputFiles = vector<string>();
  addOption("OutputFiles",&m_outputFiles, 
	    "List of the names of output files");
  
  m_varNames = vector<string>();
  addOption("VarNames",&m_varNames, 
	    "List of the names of the I/O variables");
  
  m_varScalingFactors = vector<double>();
  addOption("VarScalingFactors",&m_varScalingFactors,
            "Scaling factors for each I/O variable");
  
  m_dummyData = false;
  addOption("DummyData",&m_dummyData, 
	    "Read/write dummy data just for testing purpose");
  
  m_startIter = 0;
  addOption("StartIter",&m_startIter,
            "Starting iteration for the file processing");

  m_ioRate = 1;
  addOption("IORate",&m_ioRate,
            "Rate telling how often a file is processed", 
	    true);

}
  
//--------------------------------------------------------------------------//

FileProcessing::~FileProcessing()
{
}

//--------------------------------------------------------------------------//

void FileProcessing::configure(OptionMap& cmap, const std::string& prefix)
{ 
  LogToScreen(VERBOSE, "FileProcessing::configure() => start\n");

  BaseShockFitting::configure(cmap, prefix);
  
  // const int nbShockFitting = m_ShockFittingNames.size();
  // if (nbShockFitting > 0) {
  //   m_ShockFittingList.resize(nbShockFitting);
  // }
  
  // if (ConfigFileReader::isFirstConfig()) {
  //   for (int i = 0; i < nbShockFitting; ++i) {
  //     // create the coupling tool
  //     m_ShockFittingList[i].reset
  // 	(Factory<BaseShockFitting>::getInstance().getProvider
  // 	 (m_ShockFittingNames[i])->create(m_ShockFittingNames[i]));
  //   }
  // }
  
  // for (int i = 0; i < nbShockFitting; ++i) {
  //   // configure the coupling tool 
  //   configureDeps (cmap, m_ShockFittingList[i].get());
  // } 
  
  LogToScreen(VERBOSE, "FileProcessing::configure() => end\n");
}

//--------------------------------------------------------------------------//

void FileProcessing::getSolutionField(Field* field)
{ 
  LogToScreen(VERBOSE, "FileProcessing::getSolutionField()\n");
}

//--------------------------------------------------------------------------//

void FileProcessing::setSolutionField(Field *const field)
{ 
  LogToScreen(VERBOSE, "FileProcessing::setSolutionField()\n");
}

//--------------------------------------------------------------------------//

void FileProcessing::getMeshField(Field* field)
{ 
  LogToScreen(VERBOSE, "FileProcessing::getMeshField()\n");
}

//--------------------------------------------------------------------------//

void FileProcessing::setMeshField(Field *const field)
{ 
  LogToScreen(VERBOSE, "FileProcessing::setMeshField()\n");
}

//--------------------------------------------------------------------------//

void FileProcessing::rescaleField(Field* const field)
{
  if (m_varScalingFactors.size() > 0) {
    const unsigned nbVars = field->getStride();
    if (m_flagOutputVars.size() == 0) {
      m_flagOutputVars.resize(nbVars, true);
    }
    
    validate(m_flagOutputVars.size() == nbVars, 
	     "FileProcessing::rescaleField() => m_flagOutputVars.size() != nbVars" );
    
    unsigned countActive = 0;
    for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
      if (m_flagOutputVars[iVar]) countActive++;
    }
    
    validate(m_varScalingFactors.size() == countActive,
	     "FileProcessing::rescaleField() => m_varScalingFactors.size() != countActive");
    
    print(cout, "FileProcessing::rescaleField() => Scaling factors", m_varScalingFactors, countActive);
    
    double *const array = field->getArray();
    const unsigned fieldSize = field->getSize();
    const unsigned arraySize = fieldSize*nbVars;
    for (unsigned i = 0; i < fieldSize; ++i) {
      unsigned counter = 0;
      const unsigned start = i*nbVars; 
      for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
	if (m_flagOutputVars[iVar]) {
	  assert(start+iVar < arraySize);
	  assert(counter < m_varScalingFactors.size());
	  array[start+iVar] *= m_varScalingFactors[counter++];
	}
      }
    }
  }
}
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting
