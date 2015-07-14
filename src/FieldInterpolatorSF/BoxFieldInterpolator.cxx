// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "FieldInterpolatorSF/BoxFieldInterpolator.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Field.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<BoxFieldInterpolator, FieldInterpolator> 
boxFieldInterpolatorProv("BoxFieldInterpolator"); 
 
//--------------------------------------------------------------------------//
  
BoxFieldInterpolator::BoxFieldInterpolator(const std::string& objectName) :
  FieldInterpolator(objectName),
  m_isFirstIteration(true),
  m_dim(3),
  m_nbVarsIn(0),
  m_nbVarsOut(0),
  m_pointsInBox(),
  m_perfectMatch(),
  m_countIntp(),
  m_intpStates(),
  m_dist2AndIntp(),
  m_vminIn(),
  m_vmaxIn(),
  m_vminOut(),
  m_vmaxOut(),
  m_boxLimits()
{
  m_boxMinMax = vector<double>();
  addOption("BoxMinMax", &m_boxMinMax,
            "Coordinates of a box (xa,xb,ya,yb,za,zb) fully containing the mesh region to interpolate");
  
  m_nbSubdiv = vector<unsigned>();
  addOption("NbSubdiv", &m_nbSubdiv, "Number of box subdivisions in x,y,z");
  
  m_nbInterp = 4;
  addOption("NbInterpPoints", &m_nbInterp, "Number of interpolation points to use for each node"); 
}
  
//--------------------------------------------------------------------------//

BoxFieldInterpolator::~BoxFieldInterpolator()
{
}

//--------------------------------------------------------------------------//

void BoxFieldInterpolator::configure(SConfig::OptionMap& cmap, 
				     const std::string& prefix)
{
  FieldInterpolator::configure(cmap, prefix);
}  

//--------------------------------------------------------------------------//

void BoxFieldInterpolator::setup()
{  
  validate(m_inXyzIDs.size() > 1, 
	   "BoxFieldInterpolator::setup() => m_inXyzIDs.size() < 2");
  validate(m_outXyzIDs.size() > 1, 
	   "BoxFieldInterpolator::setup() => m_outXyzIDs.size() < 2");
  
  // set the dimension of the problem
  m_dim = std::max(m_inXyzIDs.size(), m_outXyzIDs.size());
  validate(m_dim == 2 || m_dim == 3, 
	   "BoxFieldInterpolator::setup() => m_dim != 2 && m_dim != 3");
  
  buildBoundingBox();
}

//--------------------------------------------------------------------------//

void BoxFieldInterpolator::unsetup() 
{ 
} 
 
//--------------------------------------------------------------------------// 

void BoxFieldInterpolator::buildBoundingBox()
{
  validate(m_nbSubdiv.size() > 0, "BoxFieldInterpolator::buildBox() => NbSubdiv = 0");
  validate(m_boxMinMax.size() == 6, "BoxFieldInterpolator::buildBox() => BoxMinMax size != 6");
  
  // define the boxes limits (3D)
  const unsigned nbSubX = m_nbSubdiv[XX];
  const unsigned nbSubY = m_nbSubdiv[YY];
  const unsigned nbSubZ = m_nbSubdiv[ZZ];
  const unsigned nbBoxes = nbSubX*nbSubY*nbSubZ;
  const unsigned dim2 = m_dim*2;
  m_boxLimits.resize(nbBoxes*dim2);
  
  vector<double> interval(m_dim);
  for(unsigned iDim = 0; iDim < m_dim ; ++iDim) {
    interval[iDim] = (m_boxMinMax[iDim*2+1] - m_boxMinMax[iDim*2])/m_nbSubdiv[iDim];
  }
  
  for (unsigned i = 0; i < nbSubX; ++i) {
    for (unsigned j = 0; j < nbSubY; ++j) {
      for (unsigned k = 0; k < nbSubZ; ++k) {
        const unsigned globalBoxID = k*nbSubX*nbSubY + j*nbSubX + i;
        double *const ptr = &m_boxLimits[globalBoxID*dim2];
        // [x1,x2]
        boxEntry(ptr,XX,0) = m_boxMinMax[0] + i*interval[XX];
        boxEntry(ptr,XX,1) = m_boxMinMax[0] + (i+1)*interval[XX];
        // [y1,y2]
        boxEntry(ptr,YY,0) = m_boxMinMax[2] + j*interval[YY];
        boxEntry(ptr,YY,1) = m_boxMinMax[2] + (j+1)*interval[YY];
        // [z1,z2]
        boxEntry(ptr,ZZ,0) = m_boxMinMax[4] + k*interval[ZZ];
        boxEntry(ptr,ZZ,1) = m_boxMinMax[4] + (k+1)*interval[ZZ];
	
        // print(cout, "box", &m_boxLimits[globalBoxID*dim2], 6);
      }
    }
  }
  
  m_pointsInBox.resize(nbBoxes);
}
  
//--------------------------------------------------------------------------//
  
void BoxFieldInterpolator::interpolate(Field* inField, Field* outField)
{
  LogToScreen(INFO, "BoxFieldInterpolator::interpolate()\n");
  
  if (m_isFirstIteration) {
    m_nbVarsIn  = inField->getStride() - m_dim;
    validate(m_nbVarsIn > 0, "BoxFieldInterpolator::interpolate() => m_nbVarsIn == 0");
    LogToScreen(VERBOSE, "BoxFieldInterpolator::interpolate() => m_nbVarsIn = " << m_nbVarsIn << "\n");
    
    m_nbVarsOut = outField->getStride() - m_dim;
    validate(m_nbVarsOut > 0, "BoxFieldInterpolator::interpolate() => m_nbVarsOut == 0");
    LogToScreen(VERBOSE, "BoxFieldInterpolator::interpolate() => m_nbVarsOut = " << m_nbVarsOut << "\n");
    
    const int nbPoints = outField->getSize();
    LogToScreen(VERBOSE, "BoxFieldInterpolator::interpolate() => nbPoints = " << nbPoints << "\n");
    
    m_perfectMatch.resize(nbPoints,false);
    m_countIntp.resize(nbPoints,0);
    m_intpStates.resize(nbPoints*m_nbInterp*m_nbVarsIn,0.);
    m_dist2AndIntp.resize(nbPoints*m_nbInterp);
    for (unsigned i = 0; i < m_dist2AndIntp.size(); ++i) {
      m_dist2AndIntp[i].first  = 1e12;
      m_dist2AndIntp[i].second = &m_intpStates[i*m_nbVarsIn];
    }
    
    if (m_flagOutputVars.size() == 0) {
      m_flagOutputVars.resize(m_nbVarsIn, true);
    }
    validate(m_flagOutputVars.size() == m_nbVarsIn, 
	     "BoxFieldInterpolator::interpolate() => m_flagOutputVars.size() != m_nbVarsIn");
    
    // assign each grid point of the input mesh to a box
    fillBoxes(inField); 
    
    // build interpolator
    buildInterpolator(inField, outField);
    
    m_isFirstIteration = false;
  }
  finalize(inField, outField);
}
  
//--------------------------------------------------------------------------// 

void BoxFieldInterpolator::fillBoxes(const Field *const inField)
{
  LogToScreen(VERBOSE, "BoxFieldInterpolator::fillBoxes() => start\n");
  
  // I need to have the position of the solution points !!!!!
  validate(m_dim == 3, "m_dim != 3");
  const unsigned dim2 = m_dim*2;
  const unsigned stride = inField->getStride();
  const double *const oldSol = inField->getArray();
  
  // assign points of old mesh (oldXYZ) to boxes 
  const unsigned nbPoints = inField->getSize();
  const unsigned nbBoxes  = m_pointsInBox.size();
  // clean up the boxes
  for(unsigned iBox = 0; iBox < nbBoxes; ++iBox) {
    if (m_pointsInBox[iBox].size() > 0) {
      m_pointsInBox[iBox].clear();
    }
  }
  
  unsigned counter = 0;
  for (unsigned c = 0; c < m_nbVarsIn; ++c) {
    if (m_flagOutputVars[c]) {
      counter++;
    }
  }
  
  m_vminIn.resize(counter, 1e50);
  m_vmaxIn.resize(counter, -1e50);
  
  // here we make the very reasonable assumption that xyz IDs are contiguous
  const unsigned xIDin    = m_inXyzIDs[0];
  const unsigned solIDin  = xIDin + m_dim;
  for (unsigned ip = 0; ip < nbPoints; ++ip) {
    bool foundBox = false;
    // pointer to the start of the coordinates for the current point 
    const double *const coord = &oldSol[ip*stride+xIDin];
    for(unsigned iBox = 0; iBox < nbBoxes; ++iBox) {
      double *const box = &m_boxLimits[iBox*dim2];
      
      if (boxContainsPoint(coord, box)) {
        // the given state satisfies all coordinates constraints => is included in the box
	// pointer to the start of the field vector for the current point 
        m_pointsInBox[iBox].push_back(&oldSol[ip*stride]);
        foundBox = true;
	
	// check min and max in the old solution
	unsigned varIO = 0;
	const unsigned startc = ip*stride + solIDin;
	for (unsigned c = 0; c < m_nbVarsIn; ++c) {
	  if (m_flagOutputVars[c]) {
	    const unsigned idx = startc + varIO;
	    m_vminIn[varIO] = (oldSol[idx] < m_vminIn[varIO]) ? oldSol[idx] : m_vminIn[varIO]; 
	    m_vmaxIn[varIO] = (oldSol[idx] > m_vmaxIn[varIO]) ? oldSol[idx] : m_vmaxIn[varIO];
	    varIO++;
	  }
	}
	
        break;
      }
    }
    
    if (!foundBox) {
      LogToScreen(INFO, "BoxFieldInterpolator::fillBoxes() => box not found for point (" <<  
		  coord[XX] << ", " << coord[YY] << ", " << coord[ZZ] << ")\n");
    }
    validate(foundBox,"box not found => check input for BoxMinMax");
  }
   
  for (unsigned i = 0; i < m_vminIn.size(); ++i) {
    LogToScreen(INFO, "BoxFieldInterpolator::fillboxes() => v(" << i << ") = [" 
		<< m_vminIn[i] << ", " << m_vmaxIn[i] << "]\n");
  }
  
  LogToScreen(VERBOSE, "BoxFieldInterpolator::fillBoxes() => end\n");
}
  
//----------------------------------------------------------------------------//

void BoxFieldInterpolator::buildInterpolator(Field* inField, Field* outField)
{
  LogToScreen(VERBOSE, "BoxFieldInterpolator::buildInterpolator() => start\n");
  
  // only 3D
  const unsigned dim2 = m_dim*2;
  validate(m_dim == 3, "m_dim != 3");
  const unsigned nbBoxes  = m_pointsInBox.size();
  const double dist2Min = 1e-14;
  
  // nodes, states and coordinates on local grid
  const unsigned nbPoints = outField->getSize();
  const unsigned stride   = outField->getStride();
  validate(inField->getStride() >= outField->getStride(), 
	   string("BoxFieldInterpolator::buildInterpolator() => inField->getStride() [" +
		  to_str(inField->getStride())  + "] < outField->getStride() [" + 
		  to_str(outField->getStride()) + "]"));
  
  // the nb of interpolated variables from input must be <= the nb of variables in output  
  unsigned countOutVars = 0;
  for (unsigned i = 0; i < m_flagOutputVars.size(); ++i) {
    if (m_flagOutputVars[i]) countOutVars++;
  }
  validate(countOutVars <= stride, 
	   string("BoxFieldInterpolator::buildInterpolator() => countOutVars[" + 
		  to_str(countOutVars) + "] > stride[" + to_str(stride) + "]"));
  
  // here we make the very reasonably assumption that xyz IDs are contiguous
  // and that the solution follows x,y,z 
  const unsigned xIDin    = m_inXyzIDs[0];
  const unsigned xIDout   = m_outXyzIDs[0];
  const unsigned solIDin  = xIDin + m_dim;
  const unsigned solIDout = xIDout + m_dim;
  double *const sol = outField->getArray();
  
  // nodes on the local grid
  for (unsigned i = 0; i < nbPoints; ++i) {      
    // check if any of the active boxes contains the node in question
    for(unsigned iBox = 0; iBox < nbBoxes; ++iBox) {
      const unsigned nbPointsInBox = m_pointsInBox[iBox].size();
      if (nbPointsInBox > 0) {
	// this box is not empty, so it could contain the given node
	const unsigned startg = i*stride;
	double *const box   = &m_boxLimits[iBox*dim2];
	double *const coord = &sol[startg+xIDout];
	
	if (boxContainsPoint(coord, box)) {
	  // scan all the points in the box and store states corresponding 
	  // to the m_nbInterp minimum distances
	  const vector<const double*>& oldPoints = m_pointsInBox[iBox];  
	  
	  // oldPoints[p][m_inXyzIDs[*]] store pointers to input xyz (from donor mesh)
	  for (unsigned p = 0; p < nbPointsInBox; ++p) {
	    const double dist2 = computeDistance(coord, &oldPoints[p][xIDin], m_dim);
	    
	    if (dist2 > dist2Min) {
	      // proceed only if there hasn't been found a match yet
	      if (!m_perfectMatch[i]) {
		const unsigned startDist = i*m_nbInterp;
		const unsigned lastEntry = startDist+m_nbInterp-1;
		assert(lastEntry < m_dist2AndIntp.size());
		const double *const oldSolut = &oldPoints[p][solIDin];
		
		// check the last entry in the current state set and 
		// see if the new distance is smaller
		if (dist2 < m_dist2AndIntp[lastEntry].first) {
		  m_dist2AndIntp[lastEntry].first = dist2;
		  
		  double *const statec = m_dist2AndIntp[lastEntry].second;
		  // reset the values of the interpolating states 
		  for (unsigned c = 0; c < m_nbVarsIn; ++c) {
		    statec[c] = oldSolut[c];
		  }
		  
		  // count the number of entries
		  m_countIntp[i]++;
		  
		  // sort the m_nbInterp states
		  sort(&m_dist2AndIntp[startDist], &m_dist2AndIntp[startDist]+m_nbInterp);
		}
	      }
	    }
	    else {
	      m_perfectMatch[i] = true;
	      m_countIntp[i] = 1; // there should not be other updates for node "i"
	      
	      // assign the solution corresponding to the min distance
	      const unsigned startc = i*stride + solIDout;
	      const double *const oldSolut = &oldPoints[p][solIDin];
	      unsigned varIO = 0;
	      for (unsigned c = 0; c < m_nbVarsIn; ++c) {
		if (m_flagOutputVars[c]) {
		  sol[startc + varIO] = oldSolut[c];
		  varIO++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  LogToScreen(VERBOSE, "BoxFieldInterpolator::buildInterpolator() => end\n");
}
  
//--------------------------------------------------------------------------//

void BoxFieldInterpolator::finalize(Field* inField, Field* outField) 
{
  LogToScreen(VERBOSE, "BoxFieldInterpolator::finalize() => start\n");
  
  unsigned countMatches = 0;
  unsigned counter = 0;
  for (unsigned c = 0; c < m_nbVarsIn; ++c) {
    if (m_flagOutputVars[c]) {
      counter++;
    }
  }
  
  m_vminOut.resize(counter, 1e50);
  m_vmaxOut.resize(counter, -1e50);
  
  // here we make the very reasonably assumption that xyz IDs are contiguous
  // and that the solution follows x,y,z 
  const unsigned xIDout   = m_outXyzIDs[0];
  const unsigned solIDout = xIDout + m_dim;
  const unsigned stride  = outField->getStride();
  double *const sol = outField->getArray(); // this includes x,y,z
  
  const unsigned nbPoints = outField->getSize();
  for (unsigned i = 0; i < nbPoints; ++i) { 
    // some sanity checks before completing the interpolation
    if (m_perfectMatch[i]) {
      countMatches++;
      validate(m_perfectMatch[i] && m_countIntp[i] == 1, 
	       string("BoxFieldInterpolator::finalize() => m_perfectMatch but m_countIntp != 1 for node " + to_str(i)));
    }
    
    if ((!m_perfectMatch[i]) && m_countIntp[i] < m_nbInterp) {
      LogToScreen(INFO, "BoxFieldInterpolator::finalize() => (!m_perfectMatch) && " <<  m_countIntp[i] << " < " 
		  << m_nbInterp << " for node " << i << "not satisfied!\n");
      validate (m_countIntp[i] > 0, "BoxFieldInterpolator::finalize() => m_countIntp == 0 for node " + to_str(i));
    }
    
    if (!m_perfectMatch[i]) {
      const unsigned startp = i*m_nbInterp;
      const unsigned nIntp = std::min(m_nbInterp, m_countIntp[i]);
      assert(nIntp > 0);
      double sumDist2 = 0.;
      for (unsigned ip = 0; ip < nIntp; ++ip) {
	assert(startp + ip < m_dist2AndIntp.size());
	sumDist2 += 1./m_dist2AndIntp[startp + ip].first;
      }
      
      const double invSumDist2 = 1./sumDist2;
      const unsigned startc = i*stride + solIDout;
      for (unsigned ip = 0; ip < nIntp; ++ip) {
	const double dist2 = 1./m_dist2AndIntp[startp + ip].first;
	// pointer to solution array (pointer to entry following x,y,z)
	const double *const oldSol = m_dist2AndIntp[startp + ip].second;
	unsigned varIO = 0;
	for (unsigned c = 0; c < m_nbVarsIn; ++c) {
	  if (m_flagOutputVars[c]) {
	    const unsigned idx = startc + varIO;
	    // reset solution to 0 for the first interpolating point
	    if (ip == 0) { sol[idx] = 0.;}
	    sol[idx] += dist2*oldSol[c]*invSumDist2;  
	    
	     if (varIO == 10) {
	       cout << "point [" << ip << "] => dist2 = " << dist2  << ", oldSol[c] = " << oldSol[c] << ", invSumDist2 = " << invSumDist2 << "\n";
	       cout << "dist2*oldSol[c] = " << dist2*oldSol[c] << endl; 
	       cout << "dist2*oldSol[c]*invSumDist2 = " << dist2*oldSol[c]*invSumDist2 << endl;
	       cout << "sol[idx] = " << sol[idx] << endl; 
	     }
	    varIO++;
	  }
	}
      }
      
      unsigned varIO = 0;
      for (unsigned c = 0; c < m_nbVarsIn; ++c) {
	if (m_flagOutputVars[c]) {
	  const unsigned idx = startc + varIO;
	  m_vminOut[varIO] = (sol[idx] < m_vminOut[varIO]) ? sol[idx] : m_vminOut[varIO];
	  m_vmaxOut[varIO] = (sol[idx] > m_vmaxOut[varIO]) ? sol[idx] : m_vmaxOut[varIO];
	   if (varIO == 10) {
	     cout << "m_vminIn = " << m_vminIn[varIO] << ", m_vminOut = " << m_vminOut[varIO] << "\n";
	     if (m_vminOut[varIO] < m_vminIn[varIO]) abort();
	   }
	  varIO++;
	}
      }
      
    }
  }
  
  // sanity check
  for (unsigned i = 0; i < m_vminOut.size(); ++i) {
    LogToScreen(INFO, "BoxFieldInterpolator::finalize() => v(" << i << ") = [" 
		<< m_vminOut[i] << ", " << m_vmaxOut[i] << "]\n");
    
    if (m_vminOut[i] < m_vminIn[i]) {
      LogToScreen(INFO, "BoxFieldInterpolator::finalize() => ERROR: for (" << i << ") vminOut = " 
		  << m_vminOut[i] << " < vminIn = " << m_vminIn[i] << "\n");
    } 
    if (m_vmaxOut[i] > m_vmaxIn[i]) {
      LogToScreen(INFO, "BoxFieldInterpolator::finalize() => ERROR: for (" << i << ") vmaxOut = " 
		  << m_vmaxOut[i] << " > vmaxIn = " << m_vmaxIn[i] << "\n");
    }
  }
  
  LogToScreen(VERBOSE, "BoxFieldInterpolator::finalize() => end\n");
}
  
//----------------------------------------------------------------------------//
  
} // namespace ShockFitting
