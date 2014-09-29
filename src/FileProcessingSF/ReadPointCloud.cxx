// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <limits>

#include "FileProcessingSF/ReadPointCloud.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Clist.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ReadPointCloud, FileProcessing> readPointCloudProv("ReadPointCloud"); 
  
//--------------------------------------------------------------------------//
  
ReadPointCloud::ReadPointCloud(const std::string& objectName) :
  FileProcessing(objectName), 
  m_nbVars(0),
  m_nbPoints(0),
  m_iter(0),
  m_array()
{
  m_nbDummySolVars = 0;
  addOption("NbDummySolVars", &m_nbDummySolVars, 
	    "Number of dummy solution variables to be generated"); 
}
 
//--------------------------------------------------------------------------//

ReadPointCloud::~ReadPointCloud()
{
}

//--------------------------------------------------------------------------//

void ReadPointCloud::setup()
{
  LogToScreen(VERBOSE, "ReadPointCloud::setup() => start\n");
    
  if (m_dummyData) {
    const unsigned dim = 3;
    
    ofstream file(getInputFile().c_str());
    m_nbPoints = 10000;
    file << m_nbPoints << " " << dim + m_nbDummySolVars << " ";
    
    vector<string> vnames(dim + m_nbDummySolVars);
    for (unsigned iVar = 0; iVar < dim; ++iVar) {
      vnames[iVar] = "x" + to_str(iVar);
    }
    for (unsigned iVar = 0; iVar < m_nbDummySolVars; ++iVar) {
      vnames[iVar+dim] = "v" + to_str(iVar);
    }
    for (unsigned iVar = 0; iVar < vnames.size(); ++iVar) {
      file << vnames[iVar] << " ";
    }
    file << "\n";
    
    int counter = 0;
    double x = 0.;
    double y = 0.;
    for (unsigned i = 0; i < m_nbPoints; ++i, ++counter) {
      if (counter%100 == 0) { y += 0.1; x = 0.;}
      // coordinates x, y, z 
      file << x << " " << y << " " << 0. << " ";
      // solution variables 
      for (unsigned iVar = 0; iVar < m_nbDummySolVars; ++iVar) {
	file << (iVar + 1)*10. << " ";
      }
      file << "\n";
      
      x += 0.1;
    }
    file.close();
  }
  
  // initialize the local iteration
  m_iter = getStartIter();
  
  LogToScreen(VERBOSE, "ReadPointCloud::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ReadPointCloud::unsetup() 
{ 
  LogToScreen(VERBOSE, "ReadPointCloud::unsetup()\n");
} 
 
//--------------------------------------------------------------------------// 

void ReadPointCloud::process() 
{  
  LogToScreen(INFO, "ReadPointCloud::process()\n");

  if (!fileExists(getInputFile().c_str())) {
    LogToScreen(INFO,  "ReadPointCloud::process() => ERROR: file <" 
		<< getInputFile() << "> does not exist in current directory!\n"); abort(); 
  }
  
  ifstream file(getInputFile().c_str());
  file >> m_nbPoints >> m_nbVars;
  
  if (m_dummyData) {
    validate(m_nbVars == m_nbDummySolVars + 3, 
	     "ReadPointCloud::process() => m_nbVars == m_nbDummySolVars+3");
  }
  
  vector<string> vnames(m_nbVars);
  for (unsigned iVar = 0; iVar < m_nbVars; ++iVar) {
    file >> vnames[iVar];
  }
  print(cout, "ReadPointCloud::process() => Input variables", vnames, m_nbVars);
  
  validate(m_nbPoints > 0, "ReadPointCloud::process() => m_nbPoints = 0");
  validate(m_nbVars > 0  , "ReadPointCloud::process() => m_nbVars = 0");
  
  const unsigned totSize = m_nbPoints*m_nbVars;
  if (m_array.size() > 0) {
    m_array.clear();
  }
  m_array.resize(totSize);
 
  unsigned count = 0; 
  double entry = 0.;
  for (unsigned i = 0; i < m_nbPoints; ++i) {
    for (unsigned iVar = 0; iVar < m_nbVars; ++iVar, ++count) {
      file >> entry;
      assert(count < m_array.size());
      m_array[count] = entry;
    }
  }
  
  validate(m_array.size() == m_array.capacity(), 
	   "ReadPointCloud::process() => m_array.size() != m_array.capacity()");
  
  file.close();
  
  // update the iteration number
  m_iter += getIORate();     
} 
  
//--------------------------------------------------------------------------// 

} // namespace ShockFitting
