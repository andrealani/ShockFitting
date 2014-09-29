// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockFittingObj_hh
#define ShockFitting_ShockFittingObj_hh

//--------------------------------------------------------------------------//

#include "SConfig/ConfigObject.hh"
#include "SConfig/Factory.hh"

#include "Framework/VariableTransformer.hh"
#include "Framework/FieldInterpolator.hh"
#include "Framework/FileProcessing.hh"

#define PAIR_TYPE(a) SConfig::StringT<SConfig::SharedPtr<a> >

//--------------------------------------------------------------------------//

namespace ShockFitting {
  
//--------------------------------------------------------------------------//

/// This class defines a ShockFittingObj, whose task is to provide an 
/// interface to the individual tools.
/// 
/// @author Andrea Lani

class ShockFittingObj : public SConfig::Counter, 
			 public SConfig::ConfigObject {
public:
  
  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<ShockFittingObj> PROVIDER;
  
  /// Constructor 
  /// @param objectName the concrete class name
  ShockFittingObj(const std::string& objectName);
  
  /// Destructor
  virtual ~ShockFittingObj();
  
  /// Set up this object before its first use
  virtual void setup();
  
  /// Unset up this object after its last use
  virtual void unsetup();
  
  /// Run the coupling tools
  virtual void process();
  
  /// get the name of the parent
  virtual std::string getParentName() const {return getName();}
  
  /// get the variable transformer
  std::vector<PAIR_TYPE(VariableTransformer)>& getVariableTransformerList()
  {
    return m_vTransformer;
  }
  
  /// get the field interpolators list
  std::vector<PAIR_TYPE(FieldInterpolator)>& getFieldInterpolatorList()
  {
    return m_fInterpolator;
  }
  
  /// get the file processing list
  std::vector<PAIR_TYPE(FileProcessing)>& getFileProcessingList()
  {
    return m_fProcessing;
  }
  
protected:
   
  /// create a list of coupling objects
  template <typename T>
  void createList(std::vector<PAIR_TYPE(T)>& obj)
  {
    for (unsigned i = 0; i < obj.size(); ++i) {
      const std::string name = obj[i].name();
      obj[i].ptr().reset(SConfig::Factory<T>::getInstance().
			 getProvider(name)->create(name));
    }
  }
  
  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);
    
protected:
  
  /// array of variable transformers
  std::vector<PAIR_TYPE(VariableTransformer)> m_vTransformer; 
  
  /// array of field interpolators
  std::vector<PAIR_TYPE(FieldInterpolator)> m_fInterpolator;
  
  /// array of file processing
  std::vector<PAIR_TYPE(FileProcessing)> m_fProcessing;
  
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ShockFittingObj_hh
