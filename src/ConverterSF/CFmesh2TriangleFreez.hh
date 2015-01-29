// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CFmesh2TriangleFreez_hh
#define ShockFitting_CFmesh2TriangleFreez_hh

//----------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include "Framework/Converter.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "VariableTransformerSF/Prim2Param.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CFmesh2TriangleFreez, whose task is to make conversion
/// format CFmesh format to a Triangle mesh generator format.
/// The grid is obtained from a structured mesh. Each quadrangular cell
/// is subdivide through the diagonal line in 2 triangular ones.
/// The diagonal lines could be inverted beetween two steps of the shock fitting
/// code leading to convergence problems if close the boundary layer. The cells
/// close to the boundary are therefore freezed.

class CFmesh2TriangleFreez : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CFmesh2TriangleFreez(const std::string& objectName);

  /// Destructor
  virtual ~CFmesh2TriangleFreez();

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
  std::string getClassName() { return "CFmesh2TriangleFreez"; }

  /// read CFmesh format file 
  void readCFmeshFmt();

  /// resize vectors of the MeshData pattern with the new values read on CFmesh file
  /// and assign the starting pointers for arrays 2D
  void resizeVectors();

  /// read the nodcod vector from file .node
  void readNodCod();

  /// compute edges
  void computeEdges();

  /// compute edge flag
  void computeEdgesFlag();

  /// write Triangle format file
  void writeTriangleFmt();

  /// assign variables used in ReadTrianleFmt to MeshData pattern
  void setMeshData();

  /// assign variables used in ReadTrianleFmt to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of vertices for each element
  unsigned* nvt;

  /// number of edges
  unsigned nedges;

  /// number of mesh points
  std::vector<unsigned>* npoin;
  
  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// number of boundary faces
  std::vector<unsigned>* nbfac;

  /// number of boundary points
  std::vector<unsigned>* nbpoin;

  /// number of holes
  std::vector<unsigned>* nhole;

  /// mesh points status
  std::vector <double>* zroeVect;

  /// mesh points coordinates
  std::vector <double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// vector characterizing boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing elements (assignable to MeshData)
  std::vector<int>* celcelVect;

  /// code characterizing mesh points
  std::vector<int>* nodcod;
 
  /// name of current file
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

  ///
  Array2D <int>* edgnod;

  /// 
  Array2D <int>* celedg;

  ///
  std::vector <int>flag;

  // vector of TRS_NAME strings
  std::vector<std::string> namebnd;

  /// dummy variables to store array dimension
  unsigned totsize;

  /// dummy variable to store array starting pointer
  unsigned start;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param;

  /// file storing log info
  FileLogManip logfile;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CFmesh2TriangleFreez_hh
