// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/FindShockSurfaceNormals2D.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

FindShockSurfaceNormals2D::FindShockSurfaceNormals2D()
{
}

//--------------------------------------------------------------------------//

FindShockSurfaceNormals2D::~FindShockSurfaceNormals2D()
{
}

//--------------------------------------------------------------------------//

void FindShockSurfaceNormals2D::computeUnitNormals()
{
  cout << "     => FindShockSurfaceNormals2D::computeUnitNormals()\n";

  double normalsMod;

  setPhysicsData();

  normals.resize(PhysicsInfo::getnbDim(),
                 PhysicsInfo::getnbShPointsMax(),
                 PhysicsInfo::getnbShMax());

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    // the unit normal vector is differently computed between the inner
    // points and the edge points of each shock
    // the edge points are in fact special points, their normals
    // are assigned through the setUnitNormalsSpecialPoints method
    if(ISHPOIN==0 || ISHPOIN==nShockPoints->at(ISH)-1) {}
    else {
     normals(0,ISHPOIN,ISH)=-((*XYSh)(1,ISHPOIN+1,ISH)-(*XYSh)(1,ISHPOIN,ISH));
     normals(1,ISHPOIN,ISH)=((*XYSh)(0,ISHPOIN+1,ISH)-(*XYSh)(0,ISHPOIN,ISH));
    }
    normalsMod = sqrt(normals(0,ISHPOIN,ISH)*normals(0,ISHPOIN,ISH)+
                      normals(1,ISHPOIN,ISH)*normals(1,ISHPOIN,ISH));   
    normals(0,ISHPOIN,ISH) /= normalsMod;
    normals(1,ISHPOIN,ISH) /= normalsMod;
   } 
  }

  computeUnitNormalsSpecialPoints();

  // write file plotting unit normal vector
  FILE* plotNorm;
  plotNorm = fopen("log/ShockSurfaceNormals.dat","w");

  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   fprintf(plotNorm ,"%s","TITLE = Shock Surface normals\n");
   fprintf(plotNorm ,"%s","VARIABLES = X Y Z(1) Z(2) NX NY\n");
   fprintf(plotNorm ,"%s","ZONE T='sampletext', F = FEPOINT, ET = TRIANGLE ");
   fprintf(plotNorm ,"%s %5u","N = ",nShockPoints->at(ISH));
   fprintf(plotNorm ,"%s %5u %s",", E = ",nShockPoints->at(ISH)-1,"\n");
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++)
     {fprintf(plotNorm ,"%32.16E %s",(*XYSh)(K,I,ISH)," ");}
    fprintf(plotNorm ,"%s","\n");
    fprintf(plotNorm ,"%s","1  1 ");
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++)
     {fprintf(plotNorm ,"%32.16E %s",normals(K,I,ISH)," ");}
    fprintf(plotNorm ,"%s","\n");
   }
   for (unsigned I=0; I<nShockPoints->at(ISH)-1; I++) {
    fprintf(plotNorm ,"%u %s %u %s %u %s",I+1," ",I+2," ",I+1,"\n");
   }
  }
  
  fclose(plotNorm);
}

//--------------------------------------------------------------------------//

void FindShockSurfaceNormals2D::computeUnitNormalsSpecialPoints()
{
  unsigned ISH, I, IP;

  for(unsigned ISPPNTS=0;ISPPNTS<(*nSpecPoints);ISPPNTS++) {
   if(typeSpecPoints->at(ISPPNTS)=="OPX"  ||
      typeSpecPoints->at(ISPPNTS)=="IPX"  ||
      typeSpecPoints->at(ISPPNTS)=="WPNRX"  ) {
    ISH = (*SHinSPPs)(0,0,ISPPNTS)-1;
    I = (*SHinSPPs)(1,0,ISPPNTS)-1;
    IP = I*(nShockPoints->at(ISH)-1);
    normals(0,IP,ISH) = 1.;
    normals(1,IP,ISH) = 0.;     
   }
   if(typeSpecPoints->at(ISPPNTS)=="OPY" ||
      typeSpecPoints->at(ISPPNTS)=="IPY" ||
      typeSpecPoints->at(ISPPNTS)=="WPNRY"  ) {
    ISH = (*SHinSPPs)(0,0,ISPPNTS)-1;
    I = (*SHinSPPs)(1,0,ISPPNTS)-1;
    IP = I*(nShockPoints->at(ISH)-1);
    normals(0,IP,ISH) = 0.;
    normals(1,IP,ISH) = 1.;
   }
   if(typeSpecPoints->at(ISPPNTS)=="TP") {
    cout << "FindShockSurfaceNormals2D:: (!) warning => unit normal computation ";
    cout << "                                           for TP not implemented\n";
    exit(1);
   }
   if(typeSpecPoints->at(ISPPNTS)=="QP") {
    cout << "FindShockSurfaceNormals2D:: (!) warning => unit normal computation ";
    cout << "                                           for QP not implemented\n";
    exit(1);
   }
   if(typeSpecPoints->at(ISPPNTS)=="RRX") {
    cout << "FindShockSurfaceNormals2D:: (!) warning => unit normal computation ";
    cout << "                                           for RRX not implemented\n";
    exit(1);
   }
   if(typeSpecPoints->at(ISPPNTS)=="EP") {
    cout << "FindShockSurfaceNormals2D:: (!) warning => unit normal computation ";
    cout << "                                           for EP not implemented\n";
    exit(1);
   }
   if(typeSpecPoints->at(ISPPNTS)=="C") {
    cout << "FindShockSurfaceNormals2D:: (!) warning => unit normal computation ";
    cout << "                                           for C not implemented\n";
    exit(1);
   } 
  }
}

//--------------------------------------------------------------------------//

void FindShockSurfaceNormals2D::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  nShockPoints = 
   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  typeSpecPoints =
     PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
  SHinSPPs =
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
