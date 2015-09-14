// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/GnoffoShockSensor.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<GnoffoShockSensor, ShockDetector>
gnoffoShockSensorProv("GnoffoShockSensor");

//--------------------------------------------------------------------------//

GnoffoShockSensor::GnoffoShockSensor(const std::string& objectName) :
  ShockDetector(objectName)
{
  m_maxPressRatio = 3;
  addOption("maxPressureRatio",&m_maxPressRatio,
             "Set the maximum value of the pressre ratio");
  m_minPressRatio = 2;
  addOption("minPressureRatio",&m_minPressRatio,
             "Set the minimum value of the pressure ratio");
}

//--------------------------------------------------------------------------//

GnoffoShockSensor::~GnoffoShockSensor()
{
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::setup()
{
  LogToScreen(VERBOSE, "GnoffoShockSensor::setup() => start\n");

  LogToScreen(VERBOSE, "GnoffoShockSensor::setup() => end\n");
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::unsetup()
{
  LogToScreen(VERBOSE, "GnoffoShockSensor::unsetup()\n");
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::detect()
{
  LogToScreen(INFO, "GnoffoShockSensor::detect()\n";)
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::detect(std::vector<double>& primitiveVar)
{
  LogToScreen(INFO, "GnoffoShockSensor::detect()\n";)

  // working variables
  double varMin, varMax;
  double phi;

  setMeshData();
  setPhysicsData();

  // assign starting pointers for the working array
  setAddress();

  // resize vector and array
  setSize();

  // the sensor function is computed through the pressure values
  // using varID it is possible to change the variable
  // used to compute the sensor function 
  varID = 0;

  /// the algorithm stores the shock points distribution as if only one
  /// shock is detected. The number of shock will be evaluated after 
  (*nShocks) = 1;
  unsigned ISH = 0;
  typeSh->at(ISH) = "S";

  unsigned ISHPOIN=0;
  for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {

   /// initialize varMin and varMax
   unsigned inode = (*celnod)(0,IELEM) -1;   
   varMin = primitiveVar.at(inode*(*ndof)+varID);
   varMax = primitiveVar.at(inode*(*ndof)+varID);

   // compute the min value and the max value of varID inside the element
   for(unsigned IVERT=0;IVERT<(*nvt);IVERT++) {
    inode = (*celnod)(IVERT,IELEM)-1;

    if(primitiveVar.at(inode*(*ndof)+varID)<=varMin) {
     varMin = primitiveVar.at(inode*(*ndof)+varID);
     map_varMin(0,IELEM) = varMin;
     map_varMin(1,IELEM) = inode;     
    }
    if (primitiveVar.at(inode*(*ndof)+varID)>=varMax) {
     varMax = primitiveVar.at(inode*(*ndof)+varID);
     map_varMax(0,IELEM) = varMax;
     map_varMax(1,IELEM) = inode; 
    }
   }

   phi = map_varMax(0,IELEM)/map_varMin(0,IELEM);

   if(phi>m_maxPressRatio) {

    // check in the previous evaluated elements if the detected edge
    // has been already taken into account
    bool alreadyCounted=false;
    for(unsigned JELEM=0;JELEM<IELEM;JELEM++) {
     if(map_varMin(1,IELEM)==map_varMin(1,JELEM) &&
        map_varMax(1,IELEM)==map_varMax(1,JELEM)) {
      alreadyCounted = true; break;
     }
    }

    // if the detected edge has not already taken into account,
    // store the new shock points coordinates by interpolating
    // the coordinates of the detected edge nodes
    if(!alreadyCounted) {
     (*XYSh)(0,ISHPOIN,ISH)=0.5*((*XY)(0,map_varMax(1,IELEM))+
                                (*XY)(0,map_varMin(1,IELEM)));
     (*XYSh)(1,ISHPOIN,ISH)=0.5*((*XY)(1,map_varMax(1,IELEM))+
                                (*XY)(1,map_varMin(1,IELEM)));    
     ISHPOIN++;
    }
   }
  }

  if(ISHPOIN<2) { 
   cout << "GnoffoShockSensor::error => "; 
   cout << "number of extracted shock points is less than 2\n";
   cout << "                            try another value for shockFunctMinValue\n";
   exit(1);
  }

  nShockPoints->at(ISH) = ISHPOIN;
  nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;

  // plot shock polyline
  FILE* shockPolyLine;
  shockPolyLine = fopen("log/GnoffoShockSensor.dat","w");

  fprintf(shockPolyLine,"%s","TITLE = Detected shock polyline\n");
  fprintf(shockPolyLine,"%s","VARIABLES = \"x0\" \"x1\"\n");
  fprintf(shockPolyLine,"%s","ZONE T = \"Detected ShockPolyline\"\n");
  fprintf(shockPolyLine,"%s","STRANDID=0, SOLUTIONTIME=0 ");
  fprintf(shockPolyLine,"%s %u %s","I=",nShockPoints->at(ISH),", J=1, K=1, ZONETYPE=Ordered ");
  fprintf(shockPolyLine,"%s","DATAPACKING=POINT\n");
  fprintf(shockPolyLine,"%s","DT = (SINGLE, SINGLE)\n");

  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
   fprintf(shockPolyLine,"%22.14F",(*XYSh)(0,ISHPOIN,ISH));
   fprintf(shockPolyLine,"%22.14F",(*XYSh)(1,ISHPOIN,ISH));
   fprintf(shockPolyLine,"%s","\n");
  }

  fclose(shockPolyLine);

  // de-allocate dynamic arrays
  freeArray();

  return;
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0));
  celnod = new Array2D<int> ((*nvt),nelem->at(0),
                             &celnodVect->at(0));
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::setSize()
{
  nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  XYSh->resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax());
  typeSh->resize(PhysicsInfo::getnbShMax());
  map_varMax.resize(2,nelem->at(0));
  map_varMin.resize(2,nelem->at(0));
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::freeArray()
{
  delete XY; delete celnod;
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::setMeshData()
{
  nvt = MeshData::getInstance().getData<unsigned>("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void GnoffoShockSensor::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = 
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  typeSh = PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
