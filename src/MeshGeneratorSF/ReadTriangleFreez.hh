// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReadTriangleFreez_hh
#define ShockFitting_ReadTriangleFreez_hh

//--------------------------------------------------------------------------//

#include <sstream>
#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/FileLogManip.hh"
#include "Framework/MeshGenerator.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ReadTriangleFreez, whose task is to read a Triangle Mesh 
//  Generator file and store its data.
/// There are five file types:
/// a node file if FTYPE = "node"
/// a poly file if FTYPE = "poly"
/// an ele  file if FTYPE = "ele"
/// a neigh file if FTYPE = "neigh"
/// an edge file if FTYPE = "edge"
///
/// "firstRead" flag is defined to properly store the mesh data
/// if firstRead==1 the Triangle files are read for the first time and
/// data are stored in the arrays of the background mesh
/// in the fortran version are referred to index "0" (ex: ZROE(0)
/// if firstRead==0 data are stored in the arrays of the shocked mesh
/// in order to not overwrite the old mesh ones. 
/// In the fortran version are referred to index "1"
/// (ex: ZROE(1))
/// here they are pushed back at the end of the arrays of the background mesh
/// they and are addressed by defining starting pointers.
/// The Freez version is used for the structured grids

class ReadTriangleFreez : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ReadTriangleFreez(const std::string& objectName);

  /// Destructor
  virtual ~ReadTriangleFreez();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run the mesh reading
  virtual void generate();

  /// Run the mesh reading from a given input file
  virtual void generate(std::string);

private: // helper functions

  /// return class name
  std::string getClassName() const {return "ReadTriangleFreez";}

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

  /// assign values read by ReadTriangleFreez to MeshData
  void setMeshData ();

  /// assign values read by ReadTriangleFreez to PhysicsData
  void setPhysicsData();

private: // data

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of elements vertices
  unsigned* nvt;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of edges in the mesh
  std::vector<unsigned>* nedge;

  /// number of elements in the mesh
  std::vector<unsigned>* nelem;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// number of boundary points
  std::vector<unsigned>* nbpoin;

  /// number of bodies(holes)
  std::vector<unsigned>* nhole;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh points state (assignable to MeshData)
  std::vector <double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector <double>* coorVect;

  /// mesh boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing nodes elements
  std::vector<int>* celnodVect;

  /// vector characterizing elements
  std::vector<int>* celcelVect;

  /// vector characterizing edges
  std::vector<int>* edgptrVect;

  /// mesh points coordinates (in array storing)
  Array2D <double>* coor;

  /// mesh points state (in array storing)
  Array2D <double>* zroe;

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

  /// name of the reading current file
  std::stringstream* fname;

  /// bool checking if the fname is already read from the input.case
  /// firstRead=1 not already read from the input.case
  /// firstRead=0 read already from the input.case
  unsigned* firstRead;

  /// type of the input files 
  std::vector<std::string> m_fileTypes;

  /// defines start pointers for arrays
  unsigned start;

  /// dummy variables to store size values
  unsigned totsize;
  unsigned totsize0;

  /// current number of mesh points
  unsigned m_npoin;

  /// current number of boundary faces
  unsigned m_nbfac;

  /// current number of mesh elements
  unsigned m_nelem;

  /// current number of edges
  unsigned m_nedge;

  /// current number of holes
  unsigned m_nhole;

  /// reading file
  std::ifstream file;

  /// store informations in the log file
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

#endif //ShockFitting_ReadTriangleFreez_hh
