// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_VariableTransformer_hh
#define ShockFitting_VariableTransformer_hh

//--------------------------------------------------------------------------//

#include "Framework/BaseShockFitting.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a VariableTransformer, whose task is to transform a set 
/// of logically coupled variables (states or coordinates) into another set of 
/// variables. 
/// 
/// @author Andrea Lani

class VariableTransformer : public BaseShockFitting {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<VariableTransformer> PROVIDER;
  
  /// Constructor 
  /// @param objectName the concrete class name
  VariableTransformer(const std::string& objectName);
  
  /// Destructor
  virtual ~VariableTransformer();
   
  /// Set up this object before its first use
  virtual void setup() = 0;
  
  /// Unset up this object after its last use
  virtual void unsetup() = 0;
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
  
  /// Transform one set of variable into another
  virtual void transform () = 0;

  /// Transform one set of given variable to another
  virtual void transform(std::vector<double>*, std::vector<double>*,
                         std::vector<double>*) = 0;

  /// Gets the Class name
  static std::string getClassName() {return "VariableTransformer";}

protected: // data

  /// in variables format
  std::string m_inFmt;

  /// out variables format
  std::string m_outFmt;  

  /// model used to compute the variables transformation
  std::string m_modelTransf;

  /// additional infos used to compute variables transformation
  std::string m_addInfo;
};
  
//--------------------------------------------------------------------------//
  
} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_VariableTransformer_hh
