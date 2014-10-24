// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FxMshSps_hh
#define ShockFitting_FxMshSps_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FxMshSps, whose task is to fix the mesh around
/// special points.

class FxMshSps : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  FxMshSps(const std::string& objectName);

  /// Destructor
  virtual ~FxMshSps();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// fix mesh around special points
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "FxMshSps";}

  /// create shock nodes list
  void createShockNodesList();

  /// create downstream shock edges
  void createDownShEdges();

  /// fix mesh around special points
  /// IPX, IPY, OPX, OPY, WPNRX, WPNRY
  void fixMeshForIPorOPorWPNR(unsigned ISPPNTS);

  /// fix mesh around special point RRX
  void fixMeshForRRX(unsigned ISPPNTS);

  /// set indeces characterizing the shocks
  void setShockIndeces(unsigned ,unsigned);

  /// get shock point coordinates
  void setShockPointCoor(std::string, unsigned, unsigned);

  /// check if the shock crosses a boundary point
  void checkShockBndryEdgeCrossing();

  /// Split the existing edges of the background mesh
  void splitEdges();

  /// create 2 new edges at  boundary
  void createNewEdges(unsigned, unsigned);
  void createNewEdges(unsigned, unsigned, unsigned, unsigned);

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// assign variables used in FxMshSps to MeshData pattern
  void setMeshData();

  /// assign variables used in FxMshSps to PhysicsData pattern
  void setPhysicsData();

private: // data

  /// shock indeces variables
  std::vector <unsigned> IP;
  std::vector <unsigned> ISH;

  // i-boundary edge
  int iedg1; int iedg2;

  /// i-shockpoint coordinates
  double xsh; double ysh;

  /// values returned by findBedg class
  double s1; double s2;

  /// dummy variables to store bndfac values
  int IEDGE; int I1; int I2; int IBC;

  /// boundary faces dummy index
  unsigned ibfac;

  /// dummy arrays used to create shock nodes list
  Array2D <unsigned> ISHPlistu;
  Array2D <unsigned> ISHPlistd; 

  /// space dimension
  unsigned* ndim;

  /// max number of shocks
  unsigned* nshmax;

  /// max number of shoc points for each shock
  unsigned* npshmax;

  /// number of mesh boundary faces
  unsigned* nbfac;

  /// number of shock boundary faces
  unsigned* nbfacSh;

  /// number of vertices (=3)
  unsigned* nvt;

  /// number of mesh elements
  unsigned* nelem;

  /// number of mesh points
  unsigned* npoin;

  /// number of shocks
  unsigned* r_nShocks;

  /// number of special points
  unsigned* r_nSpecPoints;

  /// number of points for each shock
  std::vector <unsigned>* r_nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* r_nShockEdges;

  /// mesh points coordinates
  std::vector <double>* coor;

  /// type of special points
  std::vector <std::string>* r_typeSpecPoints;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// mesh points coordinates (in array storing) 
  Array2D <double>* r_XY;

  /// upstream shock points coordinates
  Array3D <double>* r_XYShu;

  /// downstream shock points coordinates
  Array3D <double>* r_XYShd;

  /// array characterizing special points
  Array3D <unsigned>* r_SHinSPPs;

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_FxMshSps_hh
