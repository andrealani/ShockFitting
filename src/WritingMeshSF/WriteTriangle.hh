// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteTriangle_hh
#define ShockFitting_WriteTriangle_hh

//--------------------------------------------------------------------------//

#include <fstream>
#include "Framework/WritingMesh.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines WriteTriangle, whose task is to write Triangle Mesh
/// Generator file using data computed by the code.
/// Two vectors are defines (M02M1 and M12M0) to remove not operative shock 
/// nodes and phantom nodes.
/// They map working memory locations. 
/// Therefore it excludes nodcod points with values equal to -2, -1, -99.

class WriteTriangle : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteTriangle(const std::string& objectName);

  /// Destructor
  virtual ~WriteTriangle();

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

  /// assign variables used in WriteTriangle to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteTriangle to PhysicsData pattern
  void setPhysicsData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

private: // data

  /// space dimension
  unsigned* ndim;

  /// number of degrees of freedom
  unsigned*  ndof;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// number of element vertices (=3)
  unsigned* nvt;

  /// number of mesh points
  unsigned* npoin;

  /// number of boundary faces
  unsigned* nbfac;

  /// number of shock boundary faces
  unsigned* nbfacSh;

  /// number of additional hole points
  unsigned* naddholes;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// number of shocks
  unsigned* w_nShocks;

  /// number of Shock points for each shock
  std::vector<unsigned>* w_nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* w_nShockEdges;

  /// coordinates of addiotnal hole points
  std::vector<double>* caddholes;  

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh point status
  std::vector <double>* zroe;

  /// mesh points coordinates
  std::vector <double>* coor;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// mesh points coordinates (in array storing)
  Array2D <double>* w_XY;

  /// code characterizing shock points
  Array2D <int>* w_NodCodSh;

  /// upstream shock points status
  Array3D <double>* w_ZRoeShu;

  /// downstream shock points status
  Array3D <double>* w_ZRoeShd;

  /// shock points coordinates
  Array3D <double>* w_XYSh;

  /// shock points coordinates belonging to upstream zone
  Array3D <double>* w_XYShu;

  /// shock points coordinates belonging to downstream zone
  Array3D <double>* w_XYShd;

  /// map vector
  std::vector <unsigned> M02M1;

  /// map vector
  std::vector <unsigned> M12M0;

  /// dummy variables 
  unsigned TNPOIN; unsigned icount; 
  unsigned ICHECK; 

  /// number of mesh holes
  unsigned nHoles;

  /// dummy variables for size map vectors setting
  unsigned ilist;

  /// ofstream variables to write triangle fmt file
  std::ofstream file;

  /// name of current file
  std::vector<std::string>* fname;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteTriangle_hh

