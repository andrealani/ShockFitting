// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Triangle2CFmesh_hh
#define ShockFitting_Triangle2CFmesh_hh

//----------------------------------------------------------------------------//

#include <string>
#include "Framework/Converter.hh"
#include "MathTools/Array2D.hh"
#include "VariableTransformerSF/Param2Prim.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines Triangle2CFmesh, whose task is to make conversion
/// from Triangle mesh generator format to CFmesh files format

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

  /// read triangle mesh generator format file
  void readTriangleFmt();

  /// write CFmesh format file
  void writeCFmeshFmt();

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

  /// max number of shock
  unsigned* nshmax;

  /// max number of points for each shock
  unsigned* npshmax;

  /// max number of edges for each shock
  unsigned* neshmax;

  /// mesh points status
  std::vector <double>* zroe;

  /// mesh points coordinates
  std::vector <double>* coor;

  /// name of the current file
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

  /// dummy vector used to store boundary edges colours
  std::vector<int> ICLR;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_param2prim;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_Triangle2CFmesh_hh
