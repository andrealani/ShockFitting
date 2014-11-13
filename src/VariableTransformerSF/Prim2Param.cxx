// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VariableTransformerSF/Prim2Param.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Prim2Param::Prim2Param(const std::string& objectName) :
 VariableTransformer(objectName)
{
}

//--------------------------------------------------------------------------//

Prim2Param::~Prim2Param()
{
  delete XY; delete zroe;
}

//--------------------------------------------------------------------------//

void Prim2Param::configure(OptionMap& cmap, const std::string& prefix)
{
  VariableTransformer::configure(cmap, prefix);
}

//--------------------------------------------------------------------------//

void Prim2Param::setAddress()
{
  totsize = npoin->at(0) + npoin->at(1) + 4 *
            PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();
  zroeVect->resize(PhysicsInfo::getnbDofMax() * totsize);
  coorVect->resize(PhysicsInfo::getnbDim() * totsize);
  
  start = PhysicsInfo::getnbDim() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  XY = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                              &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double> ((*ndof),
                               (npoin->at(1) + 2 *
                                PhysicsInfo::getnbShMax() *
                                PhysicsInfo::getnbShPointsMax()),
                                &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void Prim2Param::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void Prim2Param::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  mm = PhysicsData::getInstance().getData <vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  typemol =
    PhysicsData::getInstance().getData <vector<string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
