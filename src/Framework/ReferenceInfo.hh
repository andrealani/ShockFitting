// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ReferenceInfo_hh
#define ShockFitting_ReferenceInfo_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "SConfig/ConfigObject.hh"
#include "Framework/FileLogManip.hh"
#include "Framework/Log.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines ReferenceInfo, whose task is to read freestream
/// conditions and assign them to the PhysicsData pattern.
/// It computes:
/// .) @param gref   : isoentropic coefficient of the gas
/// .) @param Tref   : reference temperature
/// .) @param pref   : reference pressure
/// .) @param rhoref : reference density
/// .) @param uref   : reference speed
/// and uses them for adimensionalization of MODEL = TCneq

class ReferenceInfo :   public SConfig::Counter,
                        public SConfig::ConfigObject {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ReferenceInfo(const std::string& objectName);

  /// Destructor
  virtual ~ReferenceInfo();

  /// Set up object before its first use
  virtual void setup();

  ///Unset up object before its last use
  virtual void unsetup();

  /// Run reference values reading
  virtual void read();

  /// get the name of the parent
  virtual std::string getParentName() const {return getName();}

  /// get the reference species densities
  static std::vector<double> getrhoiref() { return m_rhor; }

  /// get the reference temperature
  static double getTref() { return m_Tref; }

  /// get the reference pressure
  static double getpref() { return m_pref; }

  /// get the reference speed
  static double geturef() { return m_uref; }

  /// get the reference speed direction
  static bool speedDirectionXaxis() { return m_speedDirection; }

  /// get the reference length
  static double getLref() { return m_Lref; }

  /// get the reference density
  static double getrhoref() { return m_rhoref; }

  /// get the isoentropic coefficient of the gas
  static double getgam() { return m_gam; }

  /// get the gas constant
  static double getRgas() { return m_Rgas; }

  /// set rhoref
  static void setrhoref(double rhoref) { m_rhoref = rhoref; }

protected:

  /// Configures the options for this object.
  /// To be extended by derived classes.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure(SConfig::OptionMap& cmap, const std::string& prefix);

private: // helper functions

  /// get class name
  std::string getClassName() const {return "ReferenceInfo";}

  /// set reference parameters
  void setReferenceParam();

  /// adimensionalization of TCneq
  void TCneqAdim();

  /// set @param: rhoref
  void setm_rhoref();

  /// set @params: R, alpha and Cv
  void setAlpha_Rgas_Cv();

  /// set @param: gref
  void setGamRef();

  /// set @param: gm1ref
  void setGm1Ref();

  /// assign values read by ReferenceInfo to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// isoentropic coefficient of the gas
  static double m_gam;

  /// gas constant (J/kg/K)
  static double m_Rgas;

  /// reference temperature (K)
  static double m_Tref;

  /// reference pressure (p)
  static double m_pref;

  /// freestream speed (m/s)
  static double m_uref;

  /// direction of the speed (same direction as x or not)
  static bool m_speedDirection;

  /// reference density
  static double m_rhoref;

  /// species densities
  static std::vector <double> m_rhor;

  /// reference length
  static double m_Lref;

  /// heat specific ratio (assignable to PhysicsData)
  double* gref;

  /// @param gm1ref = gref-1 (assignable to PhysicsData)
  double* gm1ref;

  /// number of species
  unsigned* nsp;

  /// species molecular weights
  std::vector <double>* mm;

  /// species heat specific ratios
  std::vector <double>* gams;

  /// names of species
  std::vector <std::string>* namesp;

  /// constant gas for each species
  std::vector <double>* Rs;

  /// formation enthalpies for each species
  std::vector <double>* hf;

  /// dummy gas constant
  double R;

  /// dummy specific heat
  double Cv;

  /// chemical species concentrations
  std::vector <double> alpha; 

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ReferenceInfo_hh
