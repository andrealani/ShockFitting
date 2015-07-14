// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/BndryFacePtrFreez.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Isortrx.hh"
#include "MathTools/Binsrc.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<BndryFacePtrFreez, Remeshing> 
 readBndryFacePtrFreezProv("BndryFacePtrFreez");

//--------------------------------------------------------------------------//

BndryFacePtrFreez::BndryFacePtrFreez(const std::string& objectName) :
 Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

BndryFacePtrFreez::~BndryFacePtrFreez()
{
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::setup()
{
  LogToScreen(VERBOSE,"BndryFacePtrFreez::setup() => start\n");

  LogToScreen(VERBOSE,"BndryFacePtrFreez::setup() => end\n");
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::unsetup()
{
  LogToScreen(VERBOSE,"BndryFacePtrFreez::unsetup()\n");
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::remesh()
{
  LogToScreen(INFO,"BndryFacePtrFreez::remesh()\n");

  setMeshData();

  setAddress();

  unsigned m_nbfac = 0;

  // define a map vector holding boundary conditions names
   vector<string> map_boundaryNames;
  if(BCmap->size()>1) {
   map_boundaryNames.resize(nbfac->at(0));
   for(unsigned IBFACE=0;IBFACE<nbfac->at(0);IBFACE++) {
    map_boundaryNames.at((*bndfac)(2,IBFACE)-1)=BCmap->at(IBFACE); 
   }
  }

  for(unsigned IEDGE=0; IEDGE<nedge->at(0); IEDGE++) {
   if( (*edgptr)(2,IEDGE) > 0 && (*edgptr)(2,IEDGE) !=999) {
    (*bndfac)(0,m_nbfac)=(*edgptr)(0,IEDGE);
    (*bndfac)(1,m_nbfac)=(*edgptr)(1,IEDGE);
    (*bndfac)(2,m_nbfac)=(*edgptr)(2,IEDGE);
    if(BCmap->size()>1) {
     BCmap->at(m_nbfac)=map_boundaryNames.at((*bndfac)(2,m_nbfac)-1);
    }
     ++m_nbfac;
   }
   else if ((*edgptr)(2,IEDGE) < 0) {
    cout << "BndryFacePtrFreez::error => negative IBC: \n";
    for(unsigned IK=0; IK<3; IK++) {
     cout << (*bndfac)(IK,m_nbfac-1) << " ";
    }
    cout << endl;
   }
  }
  
  if(m_nbfac!=nbfac->at(0)) {
   cout << "BndryFacePtrFreez::error => m_nbfac should be: " << nbfac->at(0);
   cout << "                       \nwhile it is " << m_nbfac << endl;
   exit(1);
  } 

  // de-allocate dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::setAddress()
{
  unsigned start = 0;
  unsigned totsize = nbfac->at(0) + 2 *
                                    PhysicsInfo::getnbShMax() *
                                    PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  edgptr = new Array2D<int> ((*nvt), nedge->at(0), &edgptrVect->at(start));
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::freeArray()
{
  delete bndfac; delete edgptr;
}

//--------------------------------------------------------------------------//

void BndryFacePtrFreez::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nedge = MeshData::getInstance().getData <vector<unsigned> > ("NEDGE");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  edgptrVect = MeshData::getInstance().getData <vector<int> >("EDGPTR");
  BCmap = MeshData::getInstance().getData <vector<string> >("BoundariesMap");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting_hh
