// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MeshInfo_hh
#define ShockFitting_MeshInfo_hh

//--------------------------------------------------------------------------//

#include "SConfig/StringManip.hh"
#include "Framework/Connectivity.hh"
#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines MeshInfo, whose task is to read additional
/// informations on the mesh


class MeshInfo : public MeshGenerator {
public:

  /// Constructor
  /// @param objectName the concrete class name
  MeshInfo(const std::string& objectName);

  /// Destructor
  virtual ~MeshInfo();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// Run mesh info reading
  virtual void generate();

private: //helper functions

  /// assign values read by MeshInfo to MeshData
  void setMeshData();

private: //data

  /// distance between two shock faces
  double m_eps;

  /// max non dimensional distance of phantom nodes
  double m_sndmin;

  /// length of the shock edges
  double m_dxcell;

  /// number of iterations before saving solution
  unsigned m_ibak;

  /// number of hole points
  unsigned m_naddholes;

  /// hole points coordinates
  std::vector <double> m_caddholes;

  /// number of processor
  unsigned m_nproc;

  /// distance between two shock faces
  /// (assignable to MeshData)
  double* eps;

  /// max non dimensional distance of phantom nodes
  /// (assignable to MeshData)
  double* sndmin;

  /// length of the shock edges
  /// (assignable to MeshData)
  double* dxcell;

  /// number of iterations before saving solution
  /// (assignable to MeshData)
  unsigned* ibak;

  /// number of hole points
  /// (assignable to MeshData)
  unsigned* naddholes;

  /// hole points coordinates
  /// (assignable to MeshData)
  std::vector <double>* caddholes;

  /// number of processor
  /// (assignable to MeshData)
  unsigned* nproc;

};

} //namespace ShockFitting

#endif //ShockFitting_MeshInfo_hh
