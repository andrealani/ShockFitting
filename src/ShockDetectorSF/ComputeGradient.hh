// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Shock_Fitting_ComputeGradient_hh
#define Shock_Fitting_ComputeGradient_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeGradient, whose task is to evaluate
/// the gradients in the grid-points as referenced in
/// "Computational Fluid Dynamics Principles and Applications", Jiri Blazek,
/// Chapter 5, section 5.3.4

/// @author Valentina De Amicis

class ComputeGradient {
public:

  /// Constructor
  ComputeGradient();

  /// Destructor
  ~ComputeGradient();

  /// return the value of the gradient
  /// @primData   pointer to the vector of the primitive variable
  ///             used to evaluate the gradient
  Array2D<double> getGrad(std::vector<double>* primData);

  /// return the value of the normalized gradient
  /// @primData   pointer to the vector of the primitive variable
  ///             used to evaluate the gradient
  Array2D<double> getNormalizedGrad(std::vector<double>* primData);

private: // helper functions

  /// assign variables used in ComputeGradient to MeshData
  void setMeshData();

  /// assign starting pointers to array
  void setAddress();

  /// de-allocate dynamic array
  void freeArray();

private: // data

  /// number of vertices fpr each mesh element
  unsigned* nvt;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// vector storing median dual cell volume
  std::vector<double>* medianDualCellArea;

  /// vector storing nodes belongs to the elements that surrounds the grid-poin
  std::vector<unsigned>* medianDualCellNodes;

  /// vector storing the center of each median dual cell
  std::vector<unsigned>*medianDualCellNode;
  
  /// vector storing pointers to the first surrounding element-node
  std::vector<unsigned>* medianDualCellPtr;

  /// vector characterizing grid nodes
  std::vector<int>* nodcod;
 
  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;
 
  /// mesh points coordinates (in array storing)
  Array2D<double>* XY;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D<int>* celnod;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeGradient_hh
