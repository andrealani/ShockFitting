// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FndPhPs_hh
#define ShockFitting_FndPhPs_hh

//--------------------------------------------------------------------------//

#include "Framework/FileLogManip.hh"
#include "Framework/Remeshing.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FndPhPs, whose task is to find the cells crossed 
/// by the shock and the phantom points.

class FndPhPs: public Remeshing {
public:

  /// Constructor
  FndPhPs(const std::string& objectName);

  /// Destructor
  virtual ~FndPhPs();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// find phantom nodes
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "FndPhPs";}

  /// check if a cell is crossed by the shock
  bool cellCrossed(unsigned ISH_index, unsigned ielemsh_index,
                   unsigned ielem_index);

  /// compute phantom nodes
  void setPhanPoints();

  /// count phantom nodes and boundary phantom nodes
  void countPhanPoints();

  /// reset index 
  void setIndex(unsigned ISH_index, unsigned ielemsh_index,
                  unsigned ielem_index);

  /// de-allocate dynamic arrays
  void freeArray();

  /// assign variables used in FndPhPs to MeshData
  void setMeshData();

  /// assign variables used in FndPhPs to PhysicsData
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();
  
protected: // data

  /// number of vertices of a cell (=3);
  unsigned* nvt;

  /// number of shocks
  unsigned* nShocks;

  /// number of phantom points
  unsigned* nPhanPoints;

  /// number of boundary phantom points
  unsigned* nBoundPhanPoints;

  /// number of elements in the mesh
  std::vector<unsigned>* nelem;
  
  /// number of points in the mesh
  std::vector<unsigned>* npoin;

  /// number of edge in the mesh
  std::vector<unsigned>* nedge;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh points coordinates (assignable to MeshData)
  std::vector <double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// mesh points coordinates
  Array2D <double>* xy;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// dummy variable used as index
  unsigned ISH;

  /// dummy variable used as index
  unsigned ielem;

  /// dummy variable used as index
  unsigned ielemsh;

  /// cell nodes indeces
  std::vector <unsigned> n;

  /// x-coordinates of cell vertex
  std::vector <double> xc;

  /// y-coordinates of cell vertex
  std::vector <double> yc;

  /// shock points which denote straight line
  double xs1; double ys1;
 
  /// shock points which denote straight line
  double xs2; double ys2;

  /// distance between cell vertex and shock segment for each vertex
  std::vector <double> d;

  /// store informations in the log file
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_FndPhPs_hh
