// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteTriangleFreez_hh
#define ShockFitting_WriteTriangleFreez_hh

//--------------------------------------------------------------------------//

#include <stdio.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include "Framework/WritingMesh.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines WriteTriangleFreez, whose task is to write Triangle Mesh
/// Generator file using data computed by the code.
/// Two vectors are defines (M02M1 and M12M0) to remove not operative shock 
/// nodes and phantom nodes.
/// They map working memory locations. 
/// Therefore it excludes nodcod points with values equal to -2, -1, -99.
/// The Freez version should be used for the starting structured grid

class WriteTriangleFreez : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteTriangleFreez(const std::string& objectName);

  /// Destructor
  virtual ~WriteTriangleFreez();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// Write Triangle fomat file
  virtual void write();

private: // helper functions

  /// set map vector for nodcod
  void setMapVectorForNodcod();

  /// set map vector for NodCodSh
  void setMapVectorForNodcodSh();

  /// write mesh points coordinates and status on Triangle file
  void writeMeshVariables();

  /// write upstream shock points coordinates and status on Triangle file
  void writeUpstreamStatus();

  /// write downstream shock points coordinates and status on Triangle file
  void writeDownstreamStatus();

  /// set map vector for bndfac and write bndfac on poly file
  void writeBndfac();

  /// compute number of holes
  void computenbHoles();

  /// assign variables used in WriteTriangleFreez to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteTriangleFreez to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

private: // data

  /// number of degrees of freedom
  unsigned*  ndof;

  /// number of element vertices (=3)
  unsigned* nvt;

  /// number of shock boundary faces
  unsigned* nbfacSh;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh elements
  std::vector<unsigned>* nelem;
  
  /// number of mesh points
  std::vector<unsigned>* npoin;
  
  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// number of edges
  std::vector<unsigned>* nedge;

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

  /// vector characterizing edges (assignable to MeshData)
  std::vector<int>* edgptrVect;

  /// map vector
  std::vector <unsigned>* M02M1;
 
  /// map vector
  std::vector <int>* M12M0;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// edgptr(0)(i-edge) 1° endpoint of i-edge
  /// edgptr(1)(i-edge) 2° endpoint of i-edge
  /// edgptr(2)(i-edge) boundary marker of i-edge
  Array2D <int>* edgptr;

  /// mesh points coordinates (in array storing)
  Array2D <double>* XY;

  /// mesh points state (in array storing)
  Array2D <double>* Zroe;

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

  /// dummy variables 
  unsigned TNPOIN; unsigned icount; 
  unsigned ICHECK; 

  /// number of mesh holes
  unsigned nHoles;

  /// dummy variables for size map vectors setting
  unsigned ilist;

  /// variable writing the triangle fmt file
  FILE* file;

  /// name of current file
  std::stringstream* fname;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteTriangleFreez_hh

