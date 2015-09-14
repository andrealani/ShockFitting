// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockDetector_hh
#define ShockFitting_ShockDetector_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/BaseShockFitting.hh"
#include "SConfig/SharedPtr.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {


//--------------------------------------------------------------------------//

/// This class defines a ShockDetector, whose task is to find the location
/// of shock discontinuities within a CFD solution

class ShockDetector : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<ShockDetector> PROVIDER;

  /// Constructor
  /// @param objectName the concrete class n
  ShockDetector(const std::string& objectName);

  /// Destructor
  virtual ~ShockDetector();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// process detecting tools
  virtual void detect() = 0;

  ///process detecting tools with input arguments
  virtual void detect(std::vector<double>&) = 0;

  /// Gets the Class name
  static std::string getClassName() {return "ShockDetector";}

protected: // data

  /// input variables format
  std::string m_inFmt;

  /// output variables format
  std::string m_outFmt;

  /// gas model
  std::string m_modelTransf;

  /// additional infos used to compute variables transformation
  std::string m_additionalInfo;

  /// string storing the variable transformer object used for the model
  std::string m_variableTransformer;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ShockDetector_hh

