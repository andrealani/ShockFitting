// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockFileConverter_hh
#define ShockFitting_ShockFileConverter_hh

//----------------------------------------------------------------------------//

#include <fstream>
#include <vector>
#include "Framework/Converter.hh"
#include "Framework/VariableTransformer.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

/// This class defines a ShockFileConverter, whose task is to create the
/// sh00.dat input file containing the informations the shock 
/// points
/// In order to create it, a tecplot file containg the shock polyline is asked
/// The name of the tecplot file is the one specified in the input.case
/// Initializer.ShockdatFileConverter.Inputfile = "shockexample"
/// The class works only with OPY special points

class ShockFileConverter : public Converter {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ShockFileConverter(const std::string& objectName);

  /// Destructor
  virtual ~ShockFileConverter();

  /// Setup object before its first use
  virtual void setup();

  /// Unset up object after its last use
  virtual void unsetup();

  /// create the sh00.dat file
  virtual void convert();

protected: // functions

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // data

  /// input file containing the shock polyline and the downstream state
  std::string m_shInputfile;

  /// number of variables inside the shock input file (read from the input.case)
  unsigned m_nbDof;

  /// number of shocks (read from the input.case)
  unsigned m_nbShocks;

  /// number of special points (read from the input.case)
  unsigned m_nbSpecPoints;

  /// number of shock points
  unsigned nbShockPoints;

  /// type of the special point (read from the input.case)
  std::string m_typeSpecialPoint;

  /// i-shockpoint coordinates
  std::vector <double> XYSh;

  /// vectr storing primitive variables
  std::vector <double> m_prim;

  /// vector storing upstream values
  std::vector <double> ZRoeShu;

  /// vector storing downstream values
  std::vector <double> ZRoeShd;

  // ifstream variable reading the shock input file
  std::ifstream shockdat;

  /// variable writing sh00.dat file
  FILE* shfile;

  /// command object transforming variables
  PAIR_TYPE(VariableTransformer) m_prim2param; 
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

//----------------------------------------------------------------------------//

#endif // ShockFitting_ShockFileConverter_hh


