// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CFmesh2Triangle_hh
#define ShockFitting_CFmesh2Triangle_hh

//----------------------------------------------------------------------------//

#include <string>
#include "Framework/Converter.hh"
#include "MathTools/Array2D.hh"
#include "VariableTransformerSF/Prim2Param.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CFmesh2Triangle, whose task is to make conversion
/// from CFmesh format to Triangle mesh generator format

class CFmesh2Triangle : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CFmesh2Triangle(const std::string& objectName);

  /// Destructor
  virtual ~CFmesh2Triangle();

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

  /// read CFmesh format file
  void readCFmeshFmt();

  /// resize vector and array with the new values read on CFmesh file
  void resizeArrays();

  /// fill nodcod vector
  void setNodcod();

  /// count number of boundary points
  void countnbBoundaryNodes();

  /// write Triangle format file
  void writeTriangleFmt();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// assign variables used in ReadTrianleFmt to MeshData pattern
  void setMeshData();

  /// assign variables used in ReadTrianleFmt to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// space dimension
  unsigned* ndim;

  /// number of mesh points
  unsigned* npoin;

  /// number of degrees of freedom
  unsigned* ndof;

  /// max number of degrees of freedom
  unsigned* ndofmax;

  /// number of mesh elements
  unsigned* nelem;

  /// number of boundary faces
  unsigned* nbfac;

  /// number of boundary points
  unsigned* nbpoin;

  /// number of vertices for each element
  unsigned* nvt;

  /// number of holes
  unsigned* nHole;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// max number of edges for each shock
  unsigned* neshmax;

  /// mesh points status
  std::vector <double>* zroe;

  /// mesh points coordinates
  std::vector <double>* coor;

  /// code characterizing mesh points
  std::vector<int>* nodcod;

  /// name of current file
  std::vector<std::string>* fname;

  /// mesh points coordinates (in array storing)
  Array2D <double>* c_XY;

  /// mesh points status (in array storing)
  Array2D <double>* c_Zroe;

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

  // vector of TRS_NAME strings
  std::vector<std::string> namebnd;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CFmesh2Triangle_hh
