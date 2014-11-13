// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FixMshSps_hh
#define ShockFitting_FixMshSps_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FixMshSps, whose task is to fix the mesh around
/// special points.

class FixMshSps : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  FixMshSps(const std::string& objectName);

  /// Destructor
  virtual ~FixMshSps();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// fix mesh around special points
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "FixMshSps";}

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

  /// number of shock boundary faces
  unsigned* nbfacSh;

  /// number of vertices (=3)
  unsigned* nvt;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of mesh boundary faces
  std::vector<unsigned>* nbfac;

  /// number of mesh elements
  std::vector<unsigned>* nelem;
  
  /// number of mesh points
  std::vector<unsigned>* npoin;
  
  /// number of points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// mesh points coordinates (assignable to MeshData)
  std::vector <double>* coorVect;

  /// vector characterizing boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// mesh points coordinates (in array storing) 
  Array2D <double>* XY;

  /// upstream shock points coordinates
  Array3D <double>* XYShu;

  /// downstream shock points coordinates
  Array3D <double>* XYShd;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} //namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_FixMshSps_hh
