// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_TecplotFVM2Triangle_hh
#define ShockFitting_TecplotFVM2Triangle_hh

//----------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include "Framework/Converter.hh"
#include "MathTools/Array2D.hh"
#include "VariableTransformerSF/Prim2Param.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a TecplotFVM2Triangle, whose task is to make conversion
/// from Tecplot format to Triangle mesh generator format
/// the connectivity is read from .CFmesh file because COOLFluiD changes
/// the cell nodes when converting from Tecplot format and therefore the
/// cfout.plt has different cell nodes compared to the cfin.plt
/// while the LIST_STATE is read from the Tecplot file
/// It is used for the Finite Volume Method case and it differs from the
/// Tecplot2Triangle because of the .CFmesh file format read

class TecplotFVM2Triangle : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  TecplotFVM2Triangle(const std::string& objectName);

  /// Destructor
  virtual ~TecplotFVM2Triangle();

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

  /// read Tecplot format file
  void readTecplotFmt();

  /// resize vectors of the MeshData pattern with the new values read on Tecplot file
  /// and assign the starting pointers for arrays 2D
  void resizeVectors();

  /// fill nodcod vector
  void setNodcod();

  /// count number of boundary points
  void countnbBoundaryNodes();

  /// write Triangle format file
  void writeTriangleFmt();

  /// assign variables used in ReadTrianleFmt to MeshData pattern
  void setMeshData();

  /// assign variables used in ReadTrianleFmt to PhysicsData pattern
  void setPhysicsData();

  /// de-allocate the dynamic arrays
  void freeArray();

private: // data

  /// specifies if extra values are printed in addition to the
  /// commonly used degrees of freedom
  bool m_tecplotExtraValues;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of vertices for each element
  unsigned* nvt;

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

  // vector of TRS_NAME strings
  std::vector<std::string> namebnd;

  /// dummy variables to store array dimension
  unsigned totsize;

  /// dummy variable to store array starting pointer
  unsigned start;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param;

};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_TecplotFVM2Triangle_hh
