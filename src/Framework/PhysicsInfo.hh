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

  /// get the space dimensions
  static unsigned getnbDim() { return m_ndim; }

  /// get the max number of degrees of freedom
  static unsigned getnbDofMax() { return m_ndofmax; }

  /// get the max number of shocks
  static unsigned getnbShMax() { return m_nshmax; }

  /// get the max number of shock points for each shock
  static unsigned getnbShPointsMax() { return m_npshmax; }

  /// get the max number of edges for each shock
  static unsigned getnbShEdgesMax() { return m_neshmax; }

  /// get the max number of holes
  static unsigned getnbHolesMax() { return m_naddholesmax; }

  /// get the max number of special points
  static unsigned getnbSpecPointsMax() { return m_nspmax; }

  /// get the heat specific ratio
  static double getGam() {return m_gam; }

  /// get @param gm1=gam-1
  static double getGm1() { return m_gm1; }

  /// set the space dimensions
  static void setnbDim(unsigned ndim) { m_ndim = ndim; }

protected:

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: //data

  /// space dimension
  static unsigned m_ndim;

  /// heat specific ratio
  static double m_gam;

  /// max number of degree of freedom
  static unsigned m_ndofmax;

  /// max number of shocks
  static unsigned m_nshmax;

  /// max number of shock points for each shock
  static unsigned m_npshmax;

  /// max number of shock edges for each shock
  static unsigned m_neshmax;

  /// max number of special points
  static unsigned m_nspmax;

  /// max number of holes
  static unsigned m_naddholesmax;

  /// @param gm1 = gam -1
  static double m_gm1;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif //ShockFitting_PhysicsInfo_hh
