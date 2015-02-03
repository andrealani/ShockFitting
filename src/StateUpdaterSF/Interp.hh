// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Interp_hh
#define ShockFitting_Interp_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/StateUpdater.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Interp, whose task is to update values in the phantom
/// nodes of the background mesh (referred to index (0)) using values in the
/// shocked mesh (referred to index (1))

class Interp : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Interp(const std::string& objectName);

  /// Destructor
  virtual ~Interp();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// update phantom nodes values
  virtual void update();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "Interp";}

  /// assign variables used in Interp to MeshData pattern
  void setMeshData();

  /// assign variables used in Interp to PhysicsData pattern
  void setPhysicsData();

  /// set starting pointers of array 2D and 3D
  void setAddress();

  /// find the cell which the phantom node belongs to
  void finder(unsigned);

  /// return cell which the phantom node belongs to
  unsigned getCell() const { return ielem; }

  /// return boolean variable checking if the phantom node is found
  unsigned getIfound() const { return ifound; }

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of vertices
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points
  std::vector<unsigned>* nShockPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements
  std::vector<int>* celnodVect;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// map vector
  std::vector<int>* M02M1;

  /// map vector -> M02M1+NPOIN(0)
  std::vector<int>* M02M12;

  /// mesh points coordinates referring the background mesh
  Array2D <double>* XYBkg;

  /// mesh points state referring to the background mesh
  Array2D <double>* zBkg;

  /// mesh points coordinates referring the shocked mesh
  Array2D <double>* XY;

  /// mesh points state referring to the shocked mesh
  Array2D <double>* zroe;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// shock points coordinates
  Array3D<double>* XYSh;

  /// shock points coordinates belonging to upstream zone
  Array3D <double>* XYShu;

  /// shock points coordinates belonging to downstream zone
  Array3D <double>* XYShd;

  /// i-cell including the phantom node
  unsigned ielem;

  /// variable checking if the phantom node is been found
  /// @param ifound = 1 phantom node not found
  /// @param ifound = 0 phantom node found
  unsigned ifound;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Interp_hh
