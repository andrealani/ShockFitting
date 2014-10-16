// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_NormalUnitVect_hh
#define ShockFitting_NormalUnitVect_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines NormalUnitVect, whose task is compute the normal
/// unit vectors to the shocks and discontinuties
/// in shock/discontnuity points.

class NormalUnitVect: public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  NormalUnitVect(const std::string& objectName);

  /// Destructor
  virtual ~NormalUnitVect();

  /// Set up this object before its first use
  virtual void setup() {};

  /// Unset up this object before its first use
  virtual void unsetup() {};

  /// compute normal unit vector
  virtual void remesh();

private: // helper functions

  /// return gas model
  std::string getModel() const;

  /// return gas mixture
  std::string getMixture() const;

  /// return ndof
  unsigned getNbDof() const;

  /// assign variables used in NormalUnitVect to PhysicsData
  void setPhysicsData();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// gas model
  std::vector<std::string>* model;

  /// gas mixture
  std::vector<std::string>* mixture;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_NormalUnitVect_hh
