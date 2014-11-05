// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_PhysicsInfo_hh
#define ShockFitting_PhysicsInfo_hh

//--------------------------------------------------------------------------//

#include "SConfig/ConfigObject.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines PhysicsInfo, whose task is to read constant
/// values realted to the mesh and assign them to PhysicsData pattern.

class PhysicsInfo:  public SConfig::Counter,
                    public SConfig::ConfigObject {
public:

  /// Constructor
  /// @param objectName the concrete class name
  PhysicsInfo(const std::string& objectName);

  /// Destructor
  virtual ~PhysicsInfo();

  /// Set up object before its first use
  virtual void setup();

  ///Unset up object before its last use
  virtual void unsetup();

  /// Run constant values reading
  virtual void read();

  /// get the name of the parent
  virtual std::string getParentName() const {return getName();}

protected:

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper functions

  /// assign values read by PhysicsInfo to PhysicsData
  void setPhysicsData();

private: //data

  /// space dimension
  unsigned m_ndim;

  /// heat specific ratio
  double m_gam;

  /// max number of degree of freedom
  unsigned m_ndofmax;

  /// max number of shocks
  unsigned m_nshmax;

  /// max number of shock points for each shock
  unsigned m_npshmax;

  /// max number of shock edges for each shock
  unsigned m_neshmax;

  /// max number of special points
  unsigned m_nspmax;

  /// max number of holes
  unsigned m_naddholesmax;

  /// space dimension 
  ///(assignable to PhysicsData)
  unsigned* ndim;

  /// heat specific ratio
  /// assignable to PhysicsData
  double* gam;

  /// @param gm1 = gam -1
  double* gm1;

  /// max number of degree of freedom
  /// (assignable to PhysicsData)
  unsigned* ndofmax;

  /// max number of shocks
  /// (assignable to PhysicsData)
  unsigned* nshmax;

  /// max number of shock points for each shock
  /// (assignable to PhysicsData)
  unsigned* npshmax;

  /// max number of shock edges for each shock
  /// (assignable to PhysicsData)
  unsigned* neshmax;

  /// max number of special points
  /// (assignable to PhysicsData)
  unsigned* nspmax;

  /// max number of holes
  /// (assignable to PhysicsData)
  unsigned* naddholesmax;

};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_PhysicsInfo_hh
