// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoPntDispl_hh
#define Shockfitting_CoPntDispl_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include <MathTools/Array2D.hh>
#include <MathTools/Array3D.hh>

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CoPntDispl, whose task is to split the background 
/// mesh in two regions separated by a hole containing the shock. 
/// Therefore it computes displaced positions of the shock/discontinuity
/// points.

class CoPntDispl : public Remeshing {
public:


  /// Constructor
  CoPntDispl(const std::string& objectName);

  /// Destructor
  virtual ~CoPntDispl();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up object after its first use
  virtual void unsetup();

  /// Split background mesh
  virtual void remesh();


private: // helper functions

  /// return the class name
  std::string getClassName() const {return std::string("CoPntDispl");}

  /// initialize nodcodsh values to -99
  void setNodCodSh();

  /// add second layer of shock nodes
  void addSecondShockLayer();

  /// correct the displacement in the IPX OPX WPNRX 
  /// boundary special points
  void setCoorForIPXorOPXorWPNRX(unsigned ISPPNTS);

  /// correct the displacement in the IPY OPY WPNRY
  /// boundary special points 
  void setCoorForIPYorOPYorWPNRY(unsigned ISPPNTS);

  /// correct the displacement in the TP special point
  void setCoorForTP(unsigned ISPPNTS); 

  /// correct the displacement in the RRX special point
  void setCoorForRRX(unsigned ISPPNTS);

  /// correct the displacement in the QP special point
  void setCoorForQP(unsigned ISPPNTS);

  /// correct the displacement in the EP special point
  void setCoorForEP(unsigned ISPPNTS);

  /// correct the displacement in the C special point
  void setCoorForC(unsigned ISPPNTS);

  /// superimpose points for TP special points
  void superimposeDiscPoints(std::string, unsigned, unsigned,
      		             std::string, unsigned, unsigned);

  /// return upstream or downstream status vectors
  std::vector <double> getUpXVect(unsigned, unsigned);
  std::vector <double> getUpYVect(unsigned, unsigned);
  std::vector <double> getDownXVect(unsigned, unsigned);
  std::vector <double> getDownYVect(unsigned, unsigned);


  /// move discontinuity points
  void moveDiscontinuity(unsigned, unsigned);

  /// set shock characteristic indeces
  void setShockIndeces(unsigned, unsigned);

  /// assign values used in CoPntDispl to PhysicsData values
  void setPhysicsData();

  /// assign values used in CoPntDispl to PhysicsData values
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

private: // data

  /// dummy variables for the shock normal vectors
  double dx, dy;

  /// shock indeces
  std::vector<unsigned> ISH;
  std::vector<unsigned> IP;
  unsigned I;

  /// dummy variables used to compute shocks fronts coordinates
  double tx, ty, dum;

  /// vector used for the superimposition of the points
  std::vector <double> xc;
  std::vector <double> yc;
  std::vector <double> xs;
  std::vector <double> ys;

  /// dummy variables used for shocks families setting
  double f1, f3;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// type of shock
  std::vector <std::string>* typeSh;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh points status
  std::vector <double>* zroe;

  /// mesh points coordinates
  std::vector <double>* coor;

  /// code characterizing shock points
  Array2D <int>* NodCodSh;

  /// shock points coordinates
  Array3D <double>* XYSh;

  /// shock points coordinates belonging to the upstream zone
  Array3D <double>* XYShu;

  /// shock points coordinates belonging to the downstream zone
  Array3D <double>* XYShd;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// upstream state
  Array3D <double>* ZRoeShu;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// store log file infos
  FileLogManip logfile;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_CoPntDispl_hh
