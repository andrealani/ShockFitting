// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/FndPhPs.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/Remeshing.hh"
#include "MathTools/Array3D.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Ishel.hh"
#include "RemeshingSF/Rdshp.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<FndPhPs, Remeshing> findPhPsProv("FndPhPs");

//--------------------------------------------------------------------------//

FndPhPs::FndPhPs(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FndPhPs::~FndPhPs()
{
}

//--------------------------------------------------------------------------//

void FndPhPs::setup()
{
  LogToScreen(VERBOSE, "FndPhPs::setup() => start\n");

  LogToScreen(VERBOSE, "FndPhPs::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FndPhPs::unsetup()
{
  LogToScreen(VERBOSE, "FndPhPs::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FndPhPs::remesh()
{
  LogToScreen(INFO, "FndPhPs::remesh()\n");

  setMeshData();
  setPhysicsData();
  logfile.Open(getClassName());

  setAddress();
  xc.resize(3);
  yc.resize(3);
  n.resize(3);
  d.resize(3);

  unsigned cont = 0;
  for (unsigned i=0; i<nodcod->size(); i++) {
   if (nodcod->at(i)==-1 || nodcod->at(i)==-2) {cont++;}
  }

  for (unsigned iSh=0; iSh < (*nShocks); iSh++) {
   for(unsigned iElemSh=0; iElemSh < nShockEdges->at(iSh); iElemSh++) {
    for (unsigned iElem=0; iElem < nelem->at(0); iElem++) {

     if(cellCrossed(iSh,iElemSh,iElem)) {


      // for each triangle crossed by the shock element compute the distance
      // of vertices from the shock straight line.
      Rdshp dist;
      d.at(0) = dist.getRdshp(xc.at(0),yc.at(0),xs1, ys1, xs2, ys2);
      d.at(1) = dist.getRdshp(xc.at(1),yc.at(1),xs1, ys1, xs2, ys2);
      d.at(2) = dist.getRdshp(xc.at(2),yc.at(2),xs1, ys1, xs2, ys2);

      // if distance is too small the node become a phantom node
      setPhanPoints();
     } // if   
    } // iElem
   } // iElemSh
  } // iSh

  countPhanPoints();

  logfile.Close();
}

//--------------------------------------------------------------------------//

bool FndPhPs::cellCrossed(unsigned ISH_index, unsigned ielemsh_index,
                          unsigned ielem_index)
{
  setIndex(ISH_index,ielemsh_index,ielem_index); 
  for (unsigned i=0; i<(*nvt); i++) {
   n.at(i) = (*celnod)(i,ielem)-1; //c++ indeces start from 0
   xc.at(i) = (*xy)(0,n.at(i));
   yc.at(i) = (*xy)(1,n.at(i));
  }


  xs1 = (*XYSh)(0,ielemsh,ISH);
  ys1 = (*XYSh)(1,ielemsh,ISH);
  xs2 = (*XYSh)(0,ielemsh+1,ISH);
  ys2 = (*XYSh)(1,ielemsh+1,ISH);


  Ishel isCellCrossed (xc,yc,xs1,xs2,ys1,ys2);

  if (isCellCrossed.Ishel1()==0) {
   if (isCellCrossed.Ishel2()==0) {
    return true;}
  }
  return false; 
}

//--------------------------------------------------------------------------//

void FndPhPs::setPhanPoints()
{
  if(d.at(0)>0 && d.at(0)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(0))==0)
    { nodcod->at(n.at(0))=-1; }
  if(d.at(0)>0 && d.at(0)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(0))>0) 
    { nodcod->at(n.at(0))=-2; }
  if(d.at(1)>0 && d.at(1)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(1))==0)
    { nodcod->at(n.at(1))=-1; }
  if(d.at(1)>0 && d.at(1)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(1))>0) 
    { nodcod->at(n.at(1))=-2; }
  if(d.at(2)>0 && d.at(2)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(2))==0) 
    { nodcod->at(n.at(2))=-1; }
  if(d.at(2)>0 && d.at(2)<(MeshData::getInstance().getSNDMIN()) 
     && nodcod->at(n.at(2))>0)  
    { nodcod->at(n.at(2))=-2; }
}

//--------------------------------------------------------------------------//

void FndPhPs::countPhanPoints()
{
  *nPhanPoints = 0; *nBoundPhanPoints = 0;
  for (unsigned K=0; K<npoin->at(0); K++) {
  unsigned inod = K+1; // c++ indeces start from 0
   if (nodcod->at(K)==-1 || nodcod->at(K)==-2) {
     *nPhanPoints = *nPhanPoints + 1; }
   if (nodcod->at(K)==-2) { 
     *nBoundPhanPoints = *nBoundPhanPoints + 1; }
   if (nodcod->at(K)==-1) { logfile ("Node ", inod, "has become a phantom\n");}
   if (nodcod->at(K)==-2) {
    logfile ("Node ", inod, "on the boundary has become a phantom\n");}
  } 
  logfile ("Number of Phantom nodes (incl. those on the bndry)",
           *nPhanPoints, "\n");
  if (*nBoundPhanPoints > 0) {
   logfile("Uh! Oh! there are ", *nBoundPhanPoints,
           "phantom nodes on the boundary");
  }
}

//--------------------------------------------------------------------------//

void FndPhPs::setIndex(unsigned ISH_index, unsigned ielemsh_index,
                         unsigned ielem_index)
{
  ISH = ISH_index; ielem = ielem_index; ielemsh = ielemsh_index;
}

//--------------------------------------------------------------------------//

void FndPhPs::setAddress()
{
  unsigned start = 0;
  xy = new Array2D <double> (PhysicsInfo::getnbDim(),
                             npoin->at(0),
                             &coorVect->at(start));
  celnod = new Array2D <int> ((*nvt), nelem->at(0), &celnodVect->at(start));
}

//--------------------------------------------------------------------------//

void FndPhPs::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks"); 
  nShockEdges =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

void FndPhPs::setMeshData()
{
  nedge = MeshData::getInstance().getData <vector<unsigned> > ("NEDGE");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nPhanPoints = MeshData::getInstance().getData <unsigned> ("nPhanPoints");
  nBoundPhanPoints = 
             MeshData::getInstance().getData <unsigned> ("nBoundPhanPoints");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData < vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
