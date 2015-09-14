// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeShockPointsState_hh
#define ShockFitting_ComputeShockPointsState_hh

//--------------------------------------------------------------------------//

#include <algorithm>
#include <cmath>
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "Framework/FileLogManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class define a ComputeShockPointsState, whose task is to evaluate
/// the downstream and upstream state of the shock points by interpolating
/// them in the directions of the shock points normal vectors

/// @author Valentina De Amicis

class ComputeShockPointsState {
public:

  /// Constructor
  ComputeShockPointsState();

  /// Destructor
  ~ComputeShockPointsState();

  /// setup all the allocated vectors and arrays
  /// @param shockLayerThickess    distance used to extract the upstream
  ///                              and downstream coordinates 
  void setup(double shockLayerThickess);

  /// unsetup all the allocated vectors and arrays
  void unsetup();

  /// extract upstream and downstream points even if the downstream
  /// upstream state will be evaluated after
  /// @param normals   array storing normal vector to the shock points
  void extractDownstreamAndUpstreamPoints(Array3D<double> normals);

  /// evaluate the downstream and upstream state for each shock point
  /// by interpolating the state of the surrounding nodes
  void interpDownstreamAndUpstreamState();

  /// assign upstream and downstream state according to the 
  /// interpolated values
  void assignDownstreamAndUpstreamState();

private: // helper functions

   /// return class name
   std::string getClassName() const { return "ComputeShockPointsState"; }

   /// set the shock layer thickness
   void setShockLayerThickness(double shLayerThick) 
    { m_shockLayerThickness = shLayerThick; }

   /// resize vectors and array
   void setSize();

   /// assign strating pointers to array
   void setAddress();

   /// assign variables used in ComputeShockPointsState to MeshData
   void setMeshData();

   /// assign variables used in ComputeShockPointsState to PhysicsData
   void setPhysicsData();

   /// de-allocate dynamic array
   void freeArray();

private: // data

   /// distance used to extract the upstream and downstream coordinates
   double m_shockLayerThickness;

   /// number of chemical specie
   unsigned* nsp;

   /// number of vertices for each mesh element
   unsigned* nvt;

   /// number of degrees of freedom
   unsigned* ndof;

   /// number of shock points
   unsigned* nShocks;

   /// number of points for each shock
   std::vector<unsigned>* nShockPoints;

   /// number of shock edges for each shock
   std::vector<unsigned>* nShockEdges;

   /// number of mesh points
   std::vector<unsigned>* npoin;

   /// number of mesh elements
   std::vector<unsigned>* nelem;

   /// mesh points state (assigbale to MeshData)
   std::vector<double>* zroeVect;

   /// mesh points coordinates(assignable to MeshData)
   std::vector<double>* coorVect;

   /// vector characterizing nodes elements (assignable to MeshData)
   std::vector<int>* celnodVect; 

   /// mesh points state (in array storing)
   Array2D<double>* zroe;

   /// mesh points coordinates (in array storing)
   Array2D<double>* XY;

   /// celnod(0)(i-elem) 1° node of i-element
   /// celnod(1)(i-elem) 2° node of i-element
   /// celnod(2)(i-elem) 3° node of i-element
   Array2D<int>* celnod;

   /// shock points coordinates
   Array3D<double>* XYSh;

   /// upstream shock points state
   Array3D<double>* ZRoeShu;

   /// downstream shock points state
   Array3D<double>* ZRoeShd;

   /// old upstream status
   Array3D <double>* ZRoeShuOld;

   /// old downstream status
   Array3D <double>* ZRoeShdOld;

   /// working array storing upstream points coordinates
   Array3D<double> XYShu;

   /// working array storing downstream points coordinates
   Array3D<double> XYShd;

   /// file storing info
   FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ComputeShockPointsState_hh

