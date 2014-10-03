// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ReadInterpolateWrite.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/ConfigFileReader.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReadInterpolateWrite, ShockFittingObj> 
readInterpolateWriteProv("ReadInterpolateWrite");

//--------------------------------------------------------------------------//
  
ReadInterpolateWrite::ReadInterpolateWrite(const std::string& objectName) : 
  ShockFittingObj(objectName),
  m_readInputFile1(),
  m_readInputFile2(),
  m_writeOutputFile(),
  m_ioInterpolator()
{
}
  
//--------------------------------------------------------------------------//

ReadInterpolateWrite::~ReadInterpolateWrite()
{
}

//--------------------------------------------------------------------------//

void ReadInterpolateWrite::configure(SConfig::OptionMap& cmap,
				     const std::string& prefix)
{
  LogToScreen(VERBOSE, "ReadInterpolateWrite::configure() => start\n");
  
  ShockFittingObj::configure(cmap, prefix);
  
  LogToScreen(VERBOSE, "ReadInterpolateWrite::configure() => start\n");
}

//--------------------------------------------------------------------------//

void ReadInterpolateWrite::setup()
{
  LogToScreen(VERBOSE, "ReadInterpolateWrite::setup() => start\n");

  ShockFittingObj::setup();
  
  validate(m_fProcessing.size() == 3, 
	   "ReadInterpolateWrite::setup() => FileProcessingList should have size==3");
  validate(m_fInterpolator.size() == 1, 
	   "ReadInterpolateWrite::setup() => FieldInterpolatorList should have size==1");
  
  m_readInputFile1  = m_fProcessing[0].ptr();
  m_readInputFile2  = m_fProcessing[1].ptr();
  m_writeOutputFile = m_fProcessing[2].ptr();
  m_ioInterpolator  = m_fInterpolator[0].ptr();
  
  LogToScreen(VERBOSE, "ReadInterpolateWrite::setup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ReadInterpolateWrite::unsetup()
{
  LogToScreen(VERBOSE, "ReadInterpolateWrite::unsetup() => start\n");
  
  ShockFittingObj::unsetup();
  
  LogToScreen(VERBOSE, "ReadInterpolateWrite::unsetup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ReadInterpolateWrite::process()
{
  LogToScreen(VERBOSE, "ReadInterpolateWrite::process() => start\n");
  
  m_readInputFile1->process();
  m_readInputFile2->process();
  
  Field inField;
  Field outField;
  m_readInputFile1->getSolutionField(&inField);
  m_readInputFile2->getSolutionField(&outField);
  
  // interpolate from an input field to an output field
  m_ioInterpolator->interpolate(&inField, &outField);
  
  m_writeOutputFile->setSolutionField(&outField);
  m_writeOutputFile->process();
  
  LogToScreen(VERBOSE, "ReadInterpolateWrite::process() => end\n");
}
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting
