// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/IdentifySpecialPoints.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

IdentifySpecialPoints::IdentifySpecialPoints()
{
}

//--------------------------------------------------------------------------//

IdentifySpecialPoints::~IdentifySpecialPoints()
{
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::identify()
{

  cout << "     => IdentifySpecialPoints::identify()\n";

  // working variables
  double m_X, m_Y;
  Array2D<double> m_bnode(2,2);
  double middleX, middleY, distance, minDistance; 
  stringstream m_type;

  // open log file
  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  // resize vectors and array
  setSize();

  // assign starting pointers to the array
  setAddress();

  // up-to-date it works with only one shock
  assert((*nShocks)==1);

  unsigned countSpecPoint=0;
  unsigned countShEdge=0;


  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   // identify the first and the last shock points of the ISH-shock
   countShEdge=0;

   for(unsigned ISPPOINT=0;ISPPOINT<nShockPoints->at(ISH);
       ISPPOINT+=nShockPoints->at(ISH)-1) {

    assert(countShEdge<2);

    // initialize the variable which helps to find the boundary edge closest
    // to the current special point
    minDistance = MeshData::getInstance().getDXCELL()*1.5;

    m_X = (*XYSh)(0,ISPPOINT,ISH);
    m_Y = (*XYSh)(1,ISPPOINT,ISH);

    logfile("----------------------------------------\n");
    logfile("Locating shock point (",m_X, ", ",m_Y,")\n\n");

    for(unsigned IBFAC=0;IBFAC<nbfac->at(0);IBFAC++) {
     // store the coordinates of the two edge points 
     // of the IBFAC-boundary edge
     for(unsigned K=0;K<2;K++) {
      for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
       m_bnode(IDIM,K)=(*XY)(IDIM,(*bndfac)(K,IBFAC)-1); 
      }
     }

     middleX = 0.5*(m_bnode(0,0)+m_bnode(0,1));
     middleY = 0.5*(m_bnode(1,0)+m_bnode(1,1));

     // compute the distance between the special point and the middle of
     // the boundary edge
     distance = sqrt(pow((m_X-middleX),2)+pow((m_Y-middleY),2));

     // if the distance is small enough, the special point
     // is assigned to the boundary edge
     if(distance<=minDistance) {
      minDistance = distance;

      logfile("\nThe shock point is distant ", distance, "\nfrom the edge ");
      logfile("with coordinates: (",m_bnode(0,0),",",m_bnode(1,0),") and (");
      logfile(m_bnode(0,1),",",m_bnode(1,1),")\n");

      // initialize the m_type string
      m_type.str(string());

      // assign the name to the special point according to the boundary map
      // vector
      if(BCmap->at(IBFAC)=="Wall") 
       { m_type << "WPNR"; }
      else if (BCmap->at(IBFAC)=="Inlet") 
       { m_type << "IP"; }
      else if (BCmap->at(IBFAC)=="Outlet") 
       { m_type << "OP"; }
      else {
       cout << "IdentifySpecialPoints:: (!) warning => " << BCmap->at(IBFAC) << "\n";
       cout << "                        does not correspond to any special point\n";
       exit(1);
      }

      if(abs(m_bnode(0,0)-m_bnode(0,1))<=1.0e-14) {
       m_type << "Y";

       // correct the coordinates of the special point
       (*XYSh)(0,ISPPOINT,ISH)=m_bnode(0,0);
       (*XYSh)(1,ISPPOINT,ISH)=0.5*(m_bnode(1,0)+m_bnode(1,1));
//       (*XYSh)(1,ISPPOINT,ISH)=m_Y;
      }
      else if (abs(m_bnode(1,0)-m_bnode(1,1))<=1.0e-14) {
       m_type << "X";
       // correct the coordinates of the special point
       (*XYSh)(1,ISPPOINT,ISH)=m_bnode(1,0);
       (*XYSh)(0,ISPPOINT,ISH)=0.5*(m_bnode(0,0)+m_bnode(0,1));
//       (*XYSh)(0,ISPPOINT,ISH)=m_X;
      }

      logfile("whose boundary condition is \"", BCmap->at(IBFAC), "\" \n");
      logfile("assigning coordinates to the special point:");
      logfile((*XYSh)(0,ISPPOINT,ISH),",",(*XYSh)(1,ISPPOINT,ISH),"\n");
      typeSpecPoints->at(countSpecPoint) = m_type.str();
      (*SHinSPPs)(0,ISH,countSpecPoint)=ISH+1;
      (*SHinSPPs)(1,ISH,countSpecPoint)=countShEdge+1;
     }
    } // for IBFAC=0;IBFAC<nbfac->at(0)
   ++countSpecPoint; ++countShEdge;
   } // for ISPPOINT=0;ISPPOINT<nShockPoints->at(ISH)
  } // for ISH=0;ISH<(*nShocks)

  (*nSpecPoints) = countSpecPoint;

  typeSpecPoints->resize((*nSpecPoints));

  logfile ("\n-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n");
  logfile ("Recap->\n\n There are ", (*nSpecPoints), " special points in the domain:\n");
  unsigned ISHPOIN=0;
  unsigned ISH=0;
  for(unsigned ISPPNT=0;ISPPNT<(*nSpecPoints);ISPPNT++) {
   logfile(" special point # ",ISPPNT+1,"\n");
   logfile(" type: ",typeSpecPoints->at(ISPPNT),"\n");
   logfile(" leading to the ",(*SHinSPPs)(1,ISH,ISPPNT), "° edge");
   logfile(" of the ",(*SHinSPPs)(0,ISH,ISPPNT), "° shock\n");

   if      ((*SHinSPPs)(1,ISH,ISPPNT)==1) {ISHPOIN=0;}
   else if ((*SHinSPPs)(1,ISH,ISPPNT)==2) {ISHPOIN=nShockPoints->at(ISH)-1;}

   logfile(" has coordinates: ",
           (*XYSh)(0,ISHPOIN,(*SHinSPPs)(0,ISH,ISPPNT)-1),", ",
           (*XYSh)(1,ISHPOIN,(*SHinSPPs)(0,ISH,ISPPNT)-1),"\n"); 
   logfile("-----------------------------------\n");
  }

  // close the log file
  logfile.Close();

  // de-allocate dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::setSize()
{
  typeSpecPoints->resize(PhysicsInfo::getnbSpecPointsMax());
  SHinSPPs->resize(2,5,PhysicsInfo::getnbSpecPointsMax());
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0));
  totsize = nbfac->at(0) + 2 * PhysicsInfo::getnbShMax() *
                               PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(0));
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::freeArray()
{
  delete XY; delete bndfac;
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector <unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector <unsigned> > ("NBFAC");
  coorVect = MeshData::getInstance().getData <vector <double> > ("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  BCmap = MeshData::getInstance().getData <vector<string> >("BoundariesMap");
}

//--------------------------------------------------------------------------//

void IdentifySpecialPoints::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = 
   PhysicsData::getInstance().getData<vector<unsigned> > ("nShockPoints");
  nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  typeSpecPoints =
     PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  SHinSPPs =
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
  XYSh = 
       PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
