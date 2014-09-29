// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ShockFittingObj.hh"
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
ObjectProvider<ShockFittingObj, ShockFittingObj> 
ShockFittingObjProv("ShockFittingObj");

//--------------------------------------------------------------------------//
  
ShockFittingObj::ShockFittingObj(const std::string& objectName) : 
  Counter(),
  ConfigObject(objectName)
{
  m_vTransformer = vector<PAIR_TYPE(VariableTransformer)>();
  addOption("VariableTransformerList",&m_vTransformer, 
	    "List of the names of variable transformers"); 
  
  m_fInterpolator = vector<PAIR_TYPE(FieldInterpolator)>();
  addOption("FieldInterpolatorList",&m_fInterpolator, 
	    "List of the names of field interpolators"); 
  
  m_fProcessing = vector<PAIR_TYPE(FileProcessing)>();
  addOption("FileProcessingList",&m_fProcessing, 
	    "List of the names of file processing"); 
}
  
//--------------------------------------------------------------------------//

ShockFittingObj::~ShockFittingObj()
{
}

//--------------------------------------------------------------------------//

void ShockFittingObj::configure(SConfig::OptionMap& cmap,
				 const std::string& prefix)
{
  LogToScreen(VERBOSE, "ShockFittingObj::configure() => start\n");
  
  ConfigObject::configure(cmap, prefix);
  
  if (ConfigFileReader::isFirstConfig()) {
    
    // create the variable transformers
    createList<VariableTransformer>(m_vTransformer);
    
    // create the field interpolators
    createList<FieldInterpolator>(m_fInterpolator);
    
    // create the file processing
    createList<FileProcessing>(m_fProcessing);
  }
  
  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    configureDeps (cmap, m_vTransformer[i].ptr().get());
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    configureDeps (cmap, m_fInterpolator[i].ptr().get());
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    configureDeps (cmap, m_fProcessing[i].ptr().get());
  }
  
  LogToScreen(VERBOSE, "ShockFittingObj::configure() => end\n");
}

//--------------------------------------------------------------------------//

void ShockFittingObj::setup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::setup() => start\n");

  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    m_vTransformer[i].ptr()->setup();
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->setup();
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->setup();
  } 
  
  LogToScreen(VERBOSE, "ShockFittingObj::setup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::unsetup()
{
  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => start\n");
  
  // configure the variable transformers
  for (unsigned i = 0; i < m_vTransformer.size(); ++i) {
    m_vTransformer[i].ptr()->unsetup();
  }
  
  // configure the field interpolators
  for (unsigned i = 0; i < m_fInterpolator.size(); ++i) {
    m_fInterpolator[i].ptr()->unsetup();
  }
  
  // configure the file processing
  for (unsigned i = 0; i < m_fProcessing.size(); ++i) {
    m_fProcessing[i].ptr()->unsetup();
  }
  
  LogToScreen(VERBOSE, "ShockFittingObj::unsetup() => end\n");
}
  
//--------------------------------------------------------------------------//

void ShockFittingObj::process()
{
  LogToScreen(VERBOSE, "ShockFittingObj::process()\n");
}

//--------------------------------------------------------------------------//


} // namespace ShockFitting
