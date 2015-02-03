// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/ComputeStateDps.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

ComputeStateDps::ComputeStateDps(const std::string& objectName) :
  StateUpdater(objectName)
{

}

//----------------------------------------------------------------------------//

ComputeStateDps::~ComputeStateDps()
{
}

//----------------------------------------------------------------------------//

void ComputeStateDps::configure(OptionMap& cmap, const std::string& prefix)
{
  StateUpdater::configure(cmap, prefix);
}

//----------------------------------------------------------------------------//


void ComputeStateDps::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZroeShu =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeShd =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));
}

//----------------------------------------------------------------------------//

void ComputeStateDps::setDiscSpeedSize()
{
  WSh->resize(PhysicsInfo::getnbDim(),
              PhysicsInfo::getnbShPointsMax(),
              PhysicsInfo::getnbShMax());
}

//----------------------------------------------------------------------------//

void ComputeStateDps::freeArray()
{
  delete ZroeShu; delete ZroeShd;
}

//----------------------------------------------------------------------------//

void ComputeStateDps::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void ComputeStateDps::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  // this values is read in ReferenceInfo
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  typeSh = PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  ZroeShuOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZroeShdOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
