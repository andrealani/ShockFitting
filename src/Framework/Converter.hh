// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Converter_hh
#define ShockFitting_Converter_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Converter, whose task is to converter a Mesh 
/// Generator file to a file readable by the CFD solver.
/// From/To
/// .) Prim         : Primitive variables [p,u,v,T]
/// .) Cons         : Conservative variables [rho,rho_u,rho_v,rho_E]
/// .) Param        : Roe Parameter vector sqrt(rho)[1,u,v,H]
/// Model
/// .) Pg           : perfect gas
/// .) Cneq         : chemical non equilibrium
/// .) TCneq        : thermochemical non equilibrium
/// AdditionalInfos
/// for the Triangle->CFmesh conversion AdditionalInfos are
/// .) Adimensional : CF output adimensional
/// .) Dimensional  : CF output dimensional

class Converter : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<Converter> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class name
  Converter(const std::string& objectName);

  /// Destructor
  virtual ~Converter();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// convert files
  virtual void convert() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "Converter";}

protected: // data

  /// variables format input
  std::string m_inFmt;

  /// variables format output
  std::string m_outFmt;

  /// model used to compute the variables transformation
  std::string m_modelTransf;

  /// additional infos used to compute variables transformation
  std::string m_additionalInfo;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Converter_hh
