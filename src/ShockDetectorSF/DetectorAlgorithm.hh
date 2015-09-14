// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_DetectorAlgorithm_hh
#define ShockFitting_DetectorAlgorithm_hh

//----------------------------------------------------------------------------//

#include "Framework/ShockDetector.hh"
#include "Framework/FileLogManip.hh"
#include "Framework/VariableTransformer.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

class DetectorAlgorithm : public ShockDetector {
public:

  /// Constructor
  /// @param objectName the concrete class name
  DetectorAlgorithm(const std::string& objectName);

  /// Destructor
  virtual ~DetectorAlgorithm();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// command object detecting discontinuities
  virtual void detect();

  /// command object detecting discontinuities through input/output args
  /// @param firstInfoVector   vector storing info used for the shock 
  ///                          detector algorithm
  virtual void detect(std::vector<double>& firstInfoVector);

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper struct

  // sort a vector of indexes
  struct indexSort {
    const std::vector<double>& m_vec;
    indexSort(const std::vector<double>& m_vec): m_vec(m_vec) {}
    bool operator()(unsigned a, unsigned b) const { return m_vec[a] < m_vec[b]; }
  };

private: // helper functions

  /// assign variables used in the DetectorAlgorithm object to MeshData
  void setMeshData();

  /// assign variables used in the DetectorAlgorithm object to PhysicsData
  void setPhysicsData();
 
  /// assign the strating pointers to the Array2D and 3D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// specify the fitting technique
  std::string m_fittingTechnique;

  /// specify the polynomial order
  unsigned m_polynomialOrder;

  /// specify the nb of segments in the x and y direction used to
  /// split the curve
  std::vector<unsigned> m_nbSegments;

  /// specify the polynomial orders that must be used for each segment of
  /// the splitted curve
  std::vector<unsigned> m_iSegPolynomialOrders;

  /// specify if the y-coordinates are smoothed after the least squares 
  bool m_smoothOption;

  /// specify the distance to extract the upstream and downstream coodinates
  double m_shockLayerThickness;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of vertices for each mesh element
  unsigned* nvt;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// code characterizinz mesh points
  std::vector<int>* nodcod;

  /// mesh points state (assignable to MeshData)
  std::vector <double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector <double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;
 
  /// mesh points state (in array storing)
  Array2D <double>* zroe;

  /// mesh points coordinates (in array storing)
  Array2D <double>* XY;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// code characterizing shock points
  Array2D<int>* nodcodSh;

  /// shock points coordinates
  Array3D<double>* XYSh;

  /// working vector storing converted primitive variables
  std::vector<double> prim;
  
  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_param2prim;

  /// command object detecting discontinuity
  PAIR_TYPE(ShockDetector) m_whichDetector;

  /// store informations in the log file
  FileLogManip logfile;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_DetectorAlgorithm_hh
