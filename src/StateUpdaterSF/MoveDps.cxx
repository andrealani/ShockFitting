// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDps.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

MoveDps::MoveDps(const std::string& objectName) :
  StateUpdater(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDps::~MoveDps()
{
}

//----------------------------------------------------------------------------//

void MoveDps::configure(OptionMap& cmap, const std::string& prefix)
{
  StateUpdater::configure(cmap, prefix);
}

//----------------------------------------------------------------------------//

void MoveDps::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() + 
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeSh = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                 PhysicsInfo::getnbShPointsMax(),
                                 PhysicsInfo::getnbShMax(),
                                 &zroeVect->at(start));
}

//----------------------------------------------------------------------------//

void MoveDps::freeArray()
{
  delete ZroeSh;
}

//----------------------------------------------------------------------------//

void MoveDps::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void MoveDps::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  typeSh =
    PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  hf = PhysicsData::getInstance().getData <vector <double> > ("HF");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH"); 
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
