
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Remeshing_hh
#define ShockFitting_Remeshing_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/BaseShockFitting.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Remeshing, whose task is to fix the mesh created
/// by a mesh generator and do a local remesh in the neighbourhood of a
/// fitted discontinuity. 
///  
/// In particular it operates:
/// .) setting of the boundary nodes pointers
/// .) identification of phantom nodes
/// .) shock points redistribution 
/// .) splitting of the background mesh creating a hole cointaining the shock
/// .) computation of normal and tangent vectors
/// .) redefinition of the mesh around special points

class Remeshing : public BaseShockFitting {
public:

  /// typedef needed by the self-registration mechanism
  typedef SConfig::Provider<Remeshing> PROVIDER;

  /// Constructor 
  /// @param objectName the concrete class name
  Remeshing(const std::string& objectName);

  /// Destructor
  virtual ~Remeshing();

  /// Set up this object before its first use
  virtual void setup() = 0;

  /// Unset up object after its last use
  virtual void unsetup() = 0;

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

  /// process remeshing tools
  virtual void remesh() = 0;

  /// Gets the Class name
  static std::string getClassName() {return "Remeshing";}

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Remeshing_hh 
