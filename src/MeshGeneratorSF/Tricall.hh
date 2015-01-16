// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Tricall_hh
#define ShockFitting_Tricall_hh

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

//--------------------------------------------------------------------------//

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <dlfcn.h>
#include "Framework/MeshGenerator.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "TriLibrary/triangle.h"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Tricall, whose task is to create and refine a mesh
/// calling triangle as a library.

class Tricall : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Tricall(const std::string& objectName);

  /// Destructor
  virtual ~Tricall();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// create and refine the mesh
  virtual void generate();

  /// create and refine the mesh
  virtual void generate(std::string);

private: // helper functions

  /// print triangle output
  void report(struct triangulateio*,        
              int, int, int, int, int, int);

  /// assign variables used in Tricall to MeshData pattern
  void setMeshData();

  /// assign variables used in Tricall to Physics pattern
  void setPhysicsData();

  /// assign starting pointers to array 2D and 3D
  void setAddress();

  /// set sizes of Triangle vectors
  void setTriSize();

  /// set map vector for nodcod
  void setMapVectorForNodcod();

  /// set map vector for nodcodsh
  void setMapVectorForNodcodSh();

  /// store mesh points coordinates and state required by the triangulation
  void storeMeshVariables();

  /// store upstream shock points coordinates and states
  void storeUpstreamStatus();

  /// store downstream shock points coordinates and states
  void storeDownstreamStatus();

  /// set map vector for bndfac and store bndfac
  void storeBndfac();

  /// compute number of holes
  void computenbHoles();

  /// make necessary initializations so that Triangle can return
  /// a triangulation in `mid'
  void midInitialization();

  /// make necessary initializations so that Triangle can return
  /// a triangulation in `out'
  void outInitialization();

  /// Free all allocated arrays, including those allocated by Triangle
  void freeTri();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of vertices
  unsigned* nvt;

  /// number of shock points
  unsigned* nShocks;

  /// number of shock boundary faces
  unsigned* nbfacSh;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// number of Shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* nShockEdges;

  /// coordinates of additional hole points
  std::vector<double>* caddholes;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh point status (assignable to MeshData)
  std::vector <double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector <double>* coorVect;

  /// vector characterizing boundary faces
  std::vector<int>* bndfacVect;

  /// vector characterizing  nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// vector characterizing mesh elements (assignable to MeshData)
  std::vector<int>* celcelVect;

  /// map vector
  std::vector<unsigned>* M02M1;

  /// map vector
  std::vector<int>* M12M0;
  
  /// vector storing boundary faces colours
  std::vector<int>* ICLR;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// celcel(0)(i-elem) 1° neighbor of i-element
  /// celcel(1)(i-elem) 2° neighbor of i-element
  /// celcel(2)(i-elem) 3° neighbor of i-element
  Array2D <int>* celcel;

  /// mesh points coordinates (in array storing)
  Array2D <double>* XY;

  /// mesh points state (in array storing)
   Array2D <double>* zroe;

  /// code characterizing shock points
  Array2D <int>* NodCodSh;

  /// upstream shock points status
  Array3D <double>* ZRoeShu;

  /// downstream shock points status
  Array3D <double>* ZRoeShd;

  /// shock points coordinates
   Array3D <double>* XYSh;

  /// shock points coordinates belonging to upstream zone
  Array3D <double>* XYShu;

  /// shock points coordinates belonging to downstream zone
  Array3D <double>* XYShd;

  /// dummy variables used to allocate arrays
  unsigned totsize; unsigned start;

  /// dummy variables
  unsigned TNPOIN; unsigned icount;
  unsigned ICHECK;

  /// number of mesh holes
  unsigned nHoles;

  /// dummy variables for size map vectors setting
  unsigned ilist;

  /// variables used in the triangulation call
  struct triangulateio in;
  struct triangulateio mid;
  struct triangulateio out;
  struct triangulateio vorout;

  /// file storing Triangle report
  FILE* triangleReport;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Tricall_hh

