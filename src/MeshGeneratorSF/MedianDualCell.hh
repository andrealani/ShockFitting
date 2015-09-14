// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_MedianDualCell_hh
#define ShockFitting_MedianDualCell_hh

//--------------------------------------------------------------------------//

#include "Framework/MeshGenerator.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a MedianDualCell, whose task is to define the 
/// medial dual cell control volume values for each mesh point
/// as referenced in
/// "Computational Fluid Dynamics Principles and Applications", Jiri Blazek,
/// Chapter 5, section 5.1.1

/// @author Valentina De Amicis

class MedianDualCell : public MeshGenerator {
public:

  /// Constructor 
  /// @param objectName the concrete class name
  MedianDualCell(const std::string& objectName);

  /// Destructor
  virtual ~MedianDualCell();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object after its last use
  virtual void unsetup();

  /// generate the medial dual cell control volume for each grid point 
  virtual void generate();

  /// Read a given input file
  virtual void generate(std::string);

public: // helper functions

  /// return the class name
  std::string getClassName() const { return "MedianDualCell"; }

  /// erase duplicate inside a given vector
  void eraseDuplicate(std::vector<unsigned>&);

  /// assign variables used inside MedianDualCell to MeshData
  void setMeshData();

  /// assign variables used inside MedianDualCell to PhysicsData
  void setPhysicsData();

  /// assign starting pointers to array
  void setAddress();

  /// de-allocate dynamic array
  void freeArray();

private: // data

  /// number of vertices in each element
  unsigned* nvt;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// vector characterizinz mesh nodes
  std::vector<int>* nodcod;

  /// vector storing element nodes (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector storing median dual cell volume
  std::vector<double>* medianDualCellArea;

  /// vector storing center of each median dual cell
  std::vector<unsigned>* medianDualCellNode;

  /// vector storing nodes belongs to the elements that surrounds the grid-point
  std::vector<unsigned>* medianDualCellNodes;

  /// vector storing pointers to the first surrounding element-node
  std::vector<unsigned>* medianDualCellPtr;

  /// mesh points coordinates (in array storing)
  Array2D<double>* XY;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// stirng class info
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_MedianDualCell_hh
