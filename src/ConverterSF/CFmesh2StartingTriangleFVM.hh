// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CFmesh2StartingTriangleFVM_hh
#define ShockFitting_CFmesh2StartingTriangleFVM_hh

//----------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include "Framework/Converter.hh"
#include "Framework/VariableTransformer.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a CFmesh2StartingTriangleFVM, whose task is to make
/// conversion from CFmesh format to Triangle mesh generator format files
/// for the Finite Volume Method case.
/// The triangle files are the ones start the Shock Fitting algorithm
/// This class differs from the CFmesh2Triangle one because its arrays
/// do not belong to the MeshData and PhysicsData patterns

class CFmesh2StartingTriangleFVM  : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CFmesh2StartingTriangleFVM(const std::string& objectName);

  /// Destructor
  virtual ~CFmesh2StartingTriangleFVM();

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
  std::string getClassName() { return "CFmesh2StartingTriangleFVM"; }

  /// read CFmesh format file
  void readCFmeshFmt();

  /// write Triangle format file
  void writeTriangleFmt();

  /// resize vectors
  void resizeVectors();

  /// fill nodcod
  void setNodcod();

private: // data

  /// name of the CFmesh input file
  std::string m_meshInputfile;

  /// number of space dimension (read from CFmesh file)
  unsigned ndim;

  /// number of degrees of freedom (read from CFmesh file)
  unsigned ndof;

  /// number of vertices for each element
  unsigned nvt;

  /// number of mesh points (read from CFmesh file)
  unsigned npoin;

  /// number of mesh elements (read from CFmesh file)
  unsigned nelem;

  /// number of mesh boundary faces (read from CFmesh file)
  unsigned nbfac;

  /// number of holes
  unsigned nholes;

  /// mesh points state (read from the CFmesh file)
  std::vector <double> prim;

  /// mesh points state (after the variable transformer)
  std::vector <double> zroe;

  /// mesh points coordinates
  std::vector <double> XY;

  /// code characterizing mesh points
  std::vector <int> nodcod;

  /// mesh boundary faces 
  std::vector<int> bndfac;

  // vector of TRS_NAME strings
  std::vector<std::string> namebnd;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param;

  /// fstream variale reading input file
  std::ifstream file;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_CFmesh2StartingTriangleFVM_hh
