// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Triangle2CFmesh_hh
#define ShockFitting_Triangle2CFmesh_hh

//----------------------------------------------------------------------------//

#include <iomanip>
#include <string>
#include <sstream>
#include "Framework/BoundaryConnectivity.hh"
#include "Framework/Connectivity.hh"
#include "Framework/Converter.hh"
#include "Framework/Field.hh"
#include "MathTools/Array2D.hh"
#include "VariableTransformerSF/Param2Prim.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines Triangle2CFmesh, whose task is to make conversion
/// from Triangle mesh generator format to CFmesh files format.
/// The new values read on triangles files are stored in new arrays to not
/// overwrite the old mesh status.
/// These new values are pushed back at the end of the arrays of the old mesh
/// in the fortran version these new arrays are referred to index "1"
/// (ex: ZROE(1))
/// According to the Standard ShockFitting version ('original' or 
/// 'optimized') the reading and the writing objects are (or not) 
/// demanded to operate

class Triangle2CFmesh : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Triangle2CFmesh(const std::string& objectName);

  /// Destructor
  virtual ~Triangle2CFmesh();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// convert files
  virtual void convert();

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper functions

  /// return the class name
  std::string getClassName() { return "Triangle2CFmesh"; }

  /// read triangle mesh generator format file
  void readTriangleFmt();

  /// write CFmesh format file
  void writeCFmeshFmt();

  /// assign variables used in ReadTrianleFmt to MeshData pattern
  void setMeshData();

  /// assign variables used in ReadTrianleFmt to PhysicsData pattern
  void setPhysicsData();

  /// resize Connectivity array used to pass data to CF
  void resizeConnectivityArray();

  /// store CFmesh data to exchange data with CF without I/O
  void storeCFmeshData();

  /// store boundary and inner connectivity in array used to pass data to CF
  void storeConnectivity();

  /// store mesh points coordinates and state in array for CF
  void storeField();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of vertices
  unsigned* nvt;

  /// total number of boundary patches
  unsigned BNDS;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// mesh points status
  std::vector <double>* zroeVect;

  /// mesh points coordinates
  std::vector <double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// vector characterizing mesh elements (assignable to MeshData)
  std::vector<int>* celcelVect;

  /// vector characterizing boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector storing  number of boundary faces
  std::vector<int>* ICLR;

  /// face-node connectivity
  std::vector<int>* boundaryNodes;

  /// vector storing the boundary patches names
  std::vector<std::string>* boundaryNames;

  /// array of size @see BNDS x 2, storing consecutively
  /// (1) the number of elements (faces) in each boundary
  /// (2) a pointer to the first element in each boundary
  std::vector<int>* boundaryInfo;

  /// pointer to the start of the corresponding boundary face 
  std::vector<int>* boundaryPtr;

  /// pointer to the start of the corresponding element
  std::vector<int>* elementPtr;

  /// name of the current file
  std::stringstream* fname;

  /// mesh points coordinates (in array storing)
  Array2D <double>* XY;

  /// mesh points status (in array storing)
  Array2D <double>* zroe;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// celcel(0)(i-elem) 1° neighbor of i-element
  /// celcel(1)(i-elem) 2° neighbor of i-element
  /// celcel(2)(i-elem) 3° neighbor of i-element
  Array2D <int>* celcel;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// dummy variables to store array dimension
  unsigned totsize;

  /// dummy variable to store array starting pointer
  unsigned start;

  /// string for the additional info on the shock boundary
  std::string m_boundary;

  /// BoundaryConnecitivity object storing boundary mesh connectivity
  BoundaryConnectivity* inbConnectivity;

  /// Connectivity_SF object storing mesh connectivity
  Connectivity* inConnectivity;

  /// Field object storing solution variables
  Field* inStateField;

  /// Field object storing spatial coordinates
  Field* inCoordinatesField;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_param2prim;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_Triangle2CFmesh_hh
