// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeHeatFlux2DPg_hh
#define ShockFitting_computeheatFlux_hh

//--------------------------------------------------------------------------//

#include <cmath>

#include "Framework/StateUpdater.hh"
#include "Framework/VariableTransformer.hh"
#include "MathTools/Array2D.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeHeatFlux2DPg, whose task is to 
/// compute the value of the heat flux in the wall cells
/// including only the nodes that have the node used to compute the gradients 
/// leading to the same cell

/// @author Valentina De Amicis

class ComputeHeatFlux2DPg : public StateUpdater {
public:

   /// Constructor
   /// @param objectName the concrete class name
   ComputeHeatFlux2DPg(const std::string& objectName);

  /// Destructor
  virtual ~ComputeHeatFlux2DPg();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// compute the heat flux in the points of the wall cells
  virtual void update();

private: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper functions

  /// compute dynamic viscosity in the grid node using Sutherland's Law
  /// @param absT  absolute temperature in the grid node
  double getDynViscosity(double absT);

  /// compute the conductivity in the grid node
  /// @param absT absolute temperature in the grid node
  double getConductivity(double absT);

  /// assign the variables used in ComputeHeatFlux2DPg to the MeshData pattern
  void setMeshData();

  /// assign the variables used in ComputHeatFlux to the PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointer for the array2D
  void setAddress();

  /// de-allocate dynamic arrays
   void freeArray();

private: // data 

  /// number of degrees of freedoms
  unsigned* ndof;

  /// number of vertices for each cell
  unsigned* nvt;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector storing the strings of boundary conditions
  /// corresponding to the boundary marker of each boundary edge
  std::vector<std::string>* BCmap;

  /// mesh boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// mesh points state (in array storing)
  Array2D<double>* zroe;

  /// mesh points coordinate (in array storing)
  Array2D<double>* XY;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// command object making the variable conversion. From param to prim dimensional
  PAIR_TYPE(VariableTransformer) m_paramToprimDimensional;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeHeatFlux2DPg_hh

