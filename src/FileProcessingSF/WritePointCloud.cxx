// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "FileProcessingSF/WritePointCloud.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<WritePointCloud, FileProcessing> writePointCloudProv("WritePointCloud"); 
  
//--------------------------------------------------------------------------//
  
WritePointCloud::WritePointCloud(const std::string& objectName) :
  FileProcessing(objectName), 
  m_field(NULL),
  m_iter(0)
{
  m_writeTecplot = false;
  addOption("WriteTecplot",&m_writeTecplot, 
	    "Write file in TECPLOT format");
  
  m_nijk = vector<unsigned>();
  addOption("NbPointsInIJK",&m_nijk, "Number of points in I,J,K");
}
 
//--------------------------------------------------------------------------//

WritePointCloud::~WritePointCloud()
{
}

//--------------------------------------------------------------------------//

void WritePointCloud::setup()
{
  // initialize the local iteration
  m_iter = getStartIter();
}

//--------------------------------------------------------------------------//

void WritePointCloud::unsetup() 
{ 
} 
 
//--------------------------------------------------------------------------// 

void WritePointCloud::process() 
{  
  LogToScreen(INFO, "WritePointCloud::process()\n");
  
  ofstream file(getOutputFile().c_str());
  
  validate(m_field != NULL, "WritePointCloud::process() => Field not set!");
  
  const unsigned nbPoints = m_field->getSize();
  const unsigned nbInVars = m_field->getStride();
  validate(nbPoints > 0, "WritePointCloud::process() => nbPoints = 0");
  validate(nbInVars > 0, "WritePointCloud::process() => nbInVars = 0");
  
  // write a ASCII file with the following format:
  // nbPoints nbVars nameV1 nameV2 ... nameV2
  // v11 v12 ... v1n
  // v21 v22 ... v2n
  // ...
  // vn1 vn2 ... vnn
  
  unsigned nbVars = nbInVars;
  if (m_varNames.size() > 0) {
    nbVars = m_varNames.size();
    validate(nbVars <= nbInVars, "WritePointCloud::process() => nbVars > nbInVars");
  }
  
  /// the scaling factors have to be given for all I/O variables (default are =1)
  if (m_varScalingFactors.size() == 0) {
    m_varScalingFactors.resize(nbVars, 1.);
  }
  validate(m_varScalingFactors.size() == nbVars, 
	   "WritePointCloud::process() => m_varScalingFactors.size() != nbVars");
  
  vector<string> vnames(nbVars);    
  for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
    vnames[iVar] = (m_varNames.size() > 0) ? m_varNames[iVar] : "v" + to_str(iVar);
  }
  
  if (m_writeTecplot) {
    file << "TITLE = Point cloud\n" << "VARIABLES = ";
  }
  else {
    file << nbPoints << " " << nbVars << " ";
  }
  
  for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
    file << vnames[iVar] << " ";
  }
  file << "\n";
  
  if (m_writeTecplot) {
    if (m_nijk.size() == 3) {
      file << "ZONE I=" << m_nijk[XX] << ", J=" << m_nijk[YY] << ", K=" << m_nijk[ZZ] << ", DATAPACKING=POINT\n";
    }
    else {
      file << "ZONE T=\"Points\", ZONETYPE=ORDERED, DATAPACKING=POINT\n";
    }
  }
  
  print(cout, "WritePointCloud::process() => Output variables", vnames, nbVars);
  print(cout, "WritePointCloud::process() => Scaling factors", m_varScalingFactors, nbVars);
  
  if (m_flagOutputVars.size() == 0) {
    m_flagOutputVars.resize(nbInVars, true);
  }
  
  if (m_flagOutputVars.size() > 0) {
    validate(m_flagOutputVars.size() == nbInVars, 
	     "WritePointCloud::process() => m_flagOutputVars.size() != nbInVars");
    
    unsigned counter = 0;
    for (unsigned iVar = 0; iVar < nbInVars; ++iVar) {
      if (m_flagOutputVars[iVar]) counter++;
    }
    
    validate(counter == nbVars, "WritePointCloud::process() => counter != nbVars");
  }
  
  const double *const sol = m_field->getArray();
  for (unsigned i = 0; i < nbPoints; ++i) {
    const unsigned start = i*nbInVars;
    unsigned varIO = 0;
    for (unsigned iVar = 0; iVar < nbInVars; ++iVar) {
      if (m_flagOutputVars[iVar]) {
	file << sol[start + iVar]*m_varScalingFactors[varIO] << " ";
        varIO++;
      }
    }
    file << "\n";
  }

  file.close();
    
  // update the iteration number
  m_iter += getIORate();
} 
 
//--------------------------------------------------------------------------// 

// void WritePointCloud::process() 
// {  
//   LogToScreen(INFO, "WritePointCloud::process()\n");
  
//   ofstream file(getOutputFile().c_str());
  
//   validate(m_field != NULL, "WritePointCloud::process() => Field not set!");
  
//   const unsigned nbPoints = m_field->getSize();
//   const unsigned nbInVars = m_field->getStride();
//   validate(nbPoints > 0, "WritePointCloud::process() => nbPoints = 0");
//   validate(nbInVars > 0, "WritePointCloud::process() => nbInVars = 0");
  
//   // write a ASCII file with the following format:
//   // nbPoints nbVars nameV1 nameV2 ... nameV2
//   // v11 v12 ... v1n
//   // v21 v22 ... v2n
//   // ...
//   // vn1 vn2 ... vnn
  
//   unsigned nbVars = nbInVars;
//   if (m_varNames.size() > 0) {
//     nbVars = m_varNames.size();
//     validate(nbVars <= nbInVars, "WritePointCloud::process() => nbVars > nbInVars");
//   }
  
//   /// the scaling factors have to be given for all I/O variables (default are =1)
//   if (m_varScalingFactors.size() == 0) {
//     m_varScalingFactors.resize(nbVars, 1.);
//   }
//   validate(m_varScalingFactors.size() == nbVars, 
// 	   "WritePointCloud::process() => m_varScalingFactors.size() != nbVars");
  
//   vector<string> vnames(nbVars);    
//   for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
//     vnames[iVar] = (m_varNames.size() > 0) ? m_varNames[iVar] : "v" + to_str(iVar);
//   }
  
//   if (m_writeTecplot) {
//     file << "TITLE = Point cloud\n" << "VARIABLES = ";
//   }
//   else {
//     file << nbPoints << " " << nbVars << " ";
//   }
  
//   for (unsigned iVar = 0; iVar < nbVars; ++iVar) {
//     file << vnames[iVar] << " ";
//   }
//   file << "\n";
  
//   if (m_writeTecplot) {
//     if (m_nijk.size() == 3) {
//       file << "ZONE I=" << m_nijk[XX] << ", J=" << m_nijk[YY] << ", K=" << m_nijk[ZZ] << ", DATAPACKING=POINT\n";
//     }
//     else {
//       file << "ZONE T=\"Points\", ZONETYPE=ORDERED, DATAPACKING=POINT\n";
//     }
//   }
  
//   print(cout, "WritePointCloud::process() => Output variables", vnames, nbVars);
//   print(cout, "WritePointCloud::process() => Scaling factors", m_varScalingFactors, nbVars);
  
//   if (m_flagOutputVars.size() == 0) {
//     m_flagOutputVars.resize(nbInVars, true);
//   }
  
//   if (m_flagOutputVars.size() > 0) {
//     validate(m_flagOutputVars.size() == nbInVars, 
// 	     "WritePointCloud::process() => m_flagOutputVars.size() != nbInVars");
    
//     unsigned counter = 0;
//     for (unsigned iVar = 0; iVar < nbInVars; ++iVar) {
//       if (m_flagOutputVars[iVar]) counter++;
//     }
    
//     validate(counter == nbVars, "WritePointCloud::process() => counter != nbVars");
//   }
  
//   const double *const sol = m_field->getArray();
  
//   const unsigned nx = m_nijk[XX];
//   const unsigned ny = m_nijk[YY];
//   const unsigned nz = m_nijk[ZZ];
  
//   vector< vector< vector< vector<double> > > > vars(nbVars);
//   for (unsigned v = 0; v < nbVars; ++v) {
//     vars[v].resize(nx);
//     for (unsigned i = 0; i < nx; ++i) {
//       vars[v][i].resize(ny);
//       for (unsigned j = 0; j < ny; ++j) {
// 	vars[v][i][j].resize(nz);
//       }
//     }
//   }
  
//   // store in i,j,k order
//   unsigned iPoint = 0;
//   for (unsigned i = 0; i < nx; ++i) {
//     for (unsigned j= 0; j < ny; ++j) {
//       for (unsigned k = 0; k < nz; ++k, ++iPoint) {
// 	assert(iPoint < nbPoints);
// 	const unsigned start = iPoint*nbInVars;
// 	unsigned varIO = 0;
// 	for (unsigned iVar = 0; iVar < nbInVars; ++iVar) {
// 	  if (m_flagOutputVars[iVar]) {
// 	    vars[varIO][i][j][k] = sol[start + iVar]*m_varScalingFactors[varIO];
// 	    varIO++;
// 	  }
// 	}
	
//       }
//     }
//   }
  
//   // write in k,j,i order
//   iPoint = 0;
//   for (unsigned k= 0; k < nz; ++k) {
//     for (unsigned j= 0; j < ny; ++j) {
//       for (unsigned i = 0; i < nx; ++i, ++iPoint) {
// 	assert(iPoint < nbPoints);
// 	unsigned varIO = 0;
// 	for (unsigned iVar = 0; iVar < nbInVars; ++iVar) {
// 	  if (m_flagOutputVars[iVar]) {
// 	    file << vars[varIO][i][j][k] << " ";
// 	    varIO++;
// 	  }
// 	}
// 	file << "\n";
	
//       }
//     }
//   }
  
//   file.close();
    
//   // update the iteration number
//   m_iter += getIORate();
// } 
 
//--------------------------------------------------------------------------// 

} // namespace ShockFitting
