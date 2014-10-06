// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReadTriangle_hh
#define ShockFitting_ReadTriangle_hh

//--------------------------------------------------------------------------//

#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/MeshGenerator.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ReadTriangle, whose task is to read a Triangle Mesh 
//  Generator file and store its data.
/// There are five file types:
/// a node file if FTYPE = "node"
/// a poly file if FTYPE = "poly"
/// an ele  file if FTYPE = "ele"
/// a neigh file if FTYPE = "neigh"
/// an edge file if FTYPE = "edge"

class ReadTriangle : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ReadTriangle(const std::string& objectName);

  /// Destructor
  virtual ~ReadTriangle();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run the mesh reading
  virtual void generate();

private: // helper functions

  /// get the input file name
  std::string getInputFiles() const;

  /// get the node file name  
  std::string getNodeFile() const;

  /// get the poly file name
  std::string getPolyFile() const;

  /// get the ele file name
  std::string getEleFile() const;

  /// get the neigh file name
  std::string getNeighFile() const;

  /// get the edge file name
  std::string getEdgeFile() const;

  /// read .node file
  void ReadNode();

  /// read .poly file
  void ReadPoly();

  /// read .ele file
  void ReadEle();

  /// read .neigh file
  void ReadNeigh();

  /// read .neigh file
  void ReadEdge();

  /// assign values read by ReadTriangle to MeshData
  void setMeshData ();

  /// assign values read by ReadTriangle to PhysicsData
  void setPhysicsData();

private: // data

  /// type of the input files 
  std::vector<std::string> m_fileTypes;

  /// space dimension
  unsigned* ndim;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of shock elements for each shock
  unsigned* neshmax;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of mesh points
  unsigned* npoin;

  /// number of boundary points
  unsigned* nbpoin;

  /// number of boundary faces
  unsigned* nbfac;

  /// number of bodies(holes)
  unsigned* nhole;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh points status
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

  /// celcel(0)(i-elem) 1° neighbor of i-element
  /// celcel(1)(i-elem) 2° neighbor of i-element
  /// celcel(2)(i-elem) 3° neighbor of i-element
  Array2D <int>* celcel;

  /// edgptr(0)(i-edge) 1° endpoint of i-edge
  /// edgptr(1)(i-edge) 2° endpoint of i-edge
  /// edgptr(2)(i-edge) boundary marker of i-edge
  Array2D <int>* edgptr;

  /// reading file
  std::ifstream file;

};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_ReadTriangle_hh
