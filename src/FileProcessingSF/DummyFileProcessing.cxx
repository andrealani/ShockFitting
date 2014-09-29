// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "FileProcessingSF/DummyFileProcessing.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyFileProcessing, FileProcessing> 
dummyFileProcessingProv("DummyFileProcessing"); 
  
//--------------------------------------------------------------------------//
  
DummyFileProcessing::DummyFileProcessing(const std::string& objectName) :
  FileProcessing(objectName)
{
}
 
//--------------------------------------------------------------------------//

DummyFileProcessing::~DummyFileProcessing()
{
}

//--------------------------------------------------------------------------//

void DummyFileProcessing::setup()
{
}

//--------------------------------------------------------------------------//

void DummyFileProcessing::unsetup() 
{ 
} 
 
//--------------------------------------------------------------------------// 
  
void DummyFileProcessing::process()
{
  std::cout << "DummyFileProcessing::process()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
