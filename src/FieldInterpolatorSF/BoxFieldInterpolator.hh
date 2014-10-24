// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BoxFieldInterpolator_hh
#define ShockFitting_BoxFieldInterpolator_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/FieldInterpolator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a BoxFieldInterpolator, whose task is to transfer a 
/// discretized field from a grid into another, in space and/or time.
/// 
/// @author Andrea Lani

class BoxFieldInterpolator : public FieldInterpolator {
public:
  
  /// Constructor 
  /// @param objectName the concrete class name
  BoxFieldInterpolator(const std::string& objectName);
  
  /// Destructor
  virtual ~BoxFieldInterpolator();
   
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Interpolate from one field and data structure into another
  /// @param inField    input field for the solution to be interpolated
  /// @param outField   output field resulting from the interpolation inputs
  virtual void interpolate(Field* inField, Field* outField);
    
protected:
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
  /// build the bounding box
  virtual void buildBoundingBox();
  
  // assign each grid point of the old mesh to a box
  virtual void fillBoxes(const Field *const inField);
  
  /// Build the interpolator before transferring the solution from one mesh to another
  /// @param inField    input field for the solution to be interpolated
  /// @param outField   output field resulting from the interpolation inputs
  virtual void buildInterpolator(Field* inField, Field* outField);
  
  /// Finalize the field interpolation from one mesh to another
  /// @param inField    input field for the solution to be interpolated
  /// @param outField   output field resulting from the interpolation inputs
  virtual void finalize(Field* inField, Field* outField);
  
  /// access entries in the box
  double& boxEntry(double* box, int i, int j) const
  {
    return box[i*(m_dim-1)+j];
  }
  
  /// check if the box contains the given point 
  bool boxContainsPoint(const double* coord, double *const box) const
  {
    for(unsigned d = 0; d < m_dim ; ++d) {
      if(!((coord[d] >= boxEntry(box,d,0)) && (coord[d] <=  boxEntry(box,d,1)))) 
	return false;
    }
    return true;
  }
  
  /// compute the distance between two points
  double computeDistance(const double* n1, const double* n2, int size)
  {  
    double dist = 0.;
    for (int i = 0; i < size; ++i) {
      const double diff = n1[i] - n2[i];
      dist += diff*diff;
    }
    return std::sqrt(dist);
  }
  
protected:
  
  /// first iteration
  bool m_isFirstIteration;
  
  /// dimension of the problem
  unsigned m_dim;
  
  /// Number of input variables from which to interpolate 
  unsigned m_nbVarsIn;
  
  /// Number of output variables to which interpolate
  unsigned m_nbVarsOut;
  
  /// storage of points belonging to each box
  std::vector< std::vector<const double*> > m_pointsInBox;
  
  /// array of flags to indicate perfectly matching nodes that need no interpolation
  std::vector<bool> m_perfectMatch;
  
  /// array of counts for checking that the number of neighbors for each node has reached m_Nintp
  std::vector<unsigned> m_countIntp;
  
  /// array storing all states to be used for interpolation
  std::vector<double> m_intpStates;
  
  /// array storing all distances and pointers to states
  std::vector< std::pair<double, double*> > m_dist2AndIntp;
  
  /// array storing min values for each variables in the input field
  std::vector<double> m_vminIn;
  
  /// array storing max values for each variables in the input field
  std::vector<double> m_vmaxIn;
  
  /// array storing min values for each variables in the output field
  std::vector<double> m_vminOut;
  
  /// array storing max values for each variables in the output field
  std::vector<double> m_vmaxOut;
  
  /// array of logical matrices storing min-max coordinates according to the box ID
  std::vector<double> m_boxLimits;
  
  /// coordinates of a box containing the biggest mesh
  std::vector<double> m_boxMinMax;
  
  /// Number of box subdivisions in x,y,z
  std::vector<unsigned> m_nbSubdiv;
  
  /// Number of neighbor points to use for the interpolation
  unsigned m_nbInterp;
    
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

#endif // ShockFitting_BoxFieldInterpolator_hh
