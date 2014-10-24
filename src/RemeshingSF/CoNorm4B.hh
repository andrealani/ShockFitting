// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoNorm4B_hh
#define ShockFitting_CoNorm4B_hh

//--------------------------------------------------------------------------//

#include "RemeshingSF/CoNorm.hh"
#include "Framework/FileLogManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CoNorm4B, whose task is to compute the normal unit
/// vectors to the shocks and discontinuities in shock/discontinuity points
/// for NDOF=1 && MODEL="B"

class CoNorm4B : public CoNorm {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CoNorm4B(const std::string& objectName);

  /// Destructor
  ~CoNorm4B();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// Compute normal vectors
  void remesh();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("CoNorm4B");}

  /// compute tangetial vectors
  void computeTau(unsigned, unsigned);

  /// fix normal vector for typeShock = S
  void setVShNorForStype();

  /// fix normal vector for WPNRX special point
  void setVShNorForWPNRX(unsigned);

  /// fix normal vector for C special point
  void setVShNorForC(unsigned);

  /// fix normal vector for TP special point
  void setVShNorForTP(unsigned);

  /// write tecplot file
  void writeTecPlotFile();

  /// set indeces characterizing shock
  void setShockIndeces(unsigned, unsigned);

  /// take one point forward coordinates
  void onePointForward(unsigned, unsigned);

  /// take two points forward coordinates
  void twoPointsForward(unsigned, unsigned);

  /// take one point backward coordinates
  void onePointBackward(unsigned, unsigned);

  /// take two points backward coordinates
  void twoPointsBackward(unsigned, unsigned);

  /// set variables used to compute normal vectors
  void setLm();
  void setLp();
  void setTauIp1ToZero();
  void setTauIp2ToZero();
  void setTauIm1ToZero();
  void setTauIm2ToZero();

private: // data

  /// store log file infos
  FileLogManip logfile;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_CoNorm4B
