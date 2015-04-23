// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Tecplot2StartingTriangle_hh
#define ShockFitting_Tecplot2StartingTriangle_hh

//----------------------------------------------------------------------------//

#include <cmath>
#include <string>
#include <sstream>
#include "Framework/Converter.hh"
#include "Framework/VariableTransformer.hh"
#include "MathTools/Array2D.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a Tecplot2StartingTriangle, whose task is to make
/// conversion from Tecplot format to Triangle mesh generator format files
/// The triangle files are the ones start the Shock Fitting algorithm
/// This class differs from the Tecplot2Triangle one because its arrays
/// do not belong to the MeshData and PhysicsData patterns

class Tecplot2StartingTriangle  : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Tecplot2StartingTriangle(const std::string& objectName);

  /// Destructor
  virtual ~Tecplot2StartingTriangle();

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
  std::string getClassName() { return "Tecplot2StartingTriangle"; }

  /// read Tecplot format file
  void readTecplotFmt();

  /// write Triangle format file
  void writeTriangleFmt();

  /// fill nodcod
  void setNodcod();

private: // data

  /// name of the Tecplot input file
  std::vector<std::string> m_meshInputfile;

  /// number of chemical species specified in the input.case
  unsigned m_nbSpecies;

  /// specifies if extra values are printed in addition to the
  /// commonly used degrees of freedom
  bool m_tecplotExtraValues;

  /// number of degrees of freedom
  unsigned ndof;

  /// number of vertices for each element
  unsigned nvt;

  /// number of mesh points (read from Tecplot file)
  unsigned npoin;

  /// number of mesh elements (read from Tecplot file)
  unsigned nelem;

  /// number of mesh boundary faces (read from Tecplot file)
  unsigned nbfac;

  /// number of holes
  unsigned nholes;

  /// boundary colours
  unsigned NCLR;

  /// mesh points state (read from the Tecplot file)
  std::vector <double> prim;

  /// mesh points state (after the variable transformer)
  std::vector <double> zroe;

  /// mesh points coordinates
  std::vector <double> XY;

  /// code characterizing mesh points
  std::vector <int> nodcod;

  /// mesh boundary faces 
  std::vector<int> bndfac;

  /// vector of TRS_NAME strings
  std::vector<unsigned> namebnd;

  /// vector of GEOM_ENTS
  std::vector<unsigned> nFacB;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param;

  /// fstream variale reading input file
  std::ifstream file;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_Tecplot2StartingTriangle_hh
