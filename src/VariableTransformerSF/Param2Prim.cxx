// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "VariableTransformerSF/Param2Prim.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Param2Prim::Param2Prim(const std::string& objectName) :
 VariableTransformer(objectName)
{
}

//--------------------------------------------------------------------------//

Param2Prim::~Param2Prim()
{
}

//--------------------------------------------------------------------------//

void Param2Prim::configure(OptionMap& cmap, const std::string& prefix)
{
  VariableTransformer::configure(cmap, prefix);
}

//----------------------------------------------------------------------------//

void Param2Prim::setAddress()
{
  unsigned start = 0;
  v_Zroe = new Array2D <double> ((*ndof),(*npoin),&zroe->at(start));
  v_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
}

//--------------------------------------------------------------------------//

void Param2Prim::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  zroe = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coor = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void Param2Prim::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  Rgas = PhysicsData::getInstance().getData <double> ("RgasFreeStream");
  gam = PhysicsData::getInstance().getData <double> ("GamFreeStream");
  Tref = PhysicsData::getInstance().getData <double> ("TREF");
  pref = PhysicsData::getInstance().getData <double> ("PREF");
  uref = PhysicsData::getInstance().getData <double> ("UREF");
  rhoref = PhysicsData::getInstance().getData <double> ("RHOREF");
  Lref = PhysicsData::getInstance().getData <double> ("LREF");
  nmol = PhysicsData::getInstance().getData <unsigned> ("NMOL");
  ie = PhysicsData::getInstance().getData <unsigned> ("IE");
  iev = PhysicsData::getInstance().getData <unsigned> ("IEV");
  ix = PhysicsData::getInstance().getData <unsigned> ("IX");
  iy = PhysicsData::getInstance().getData <unsigned> ("IY");
  mm = PhysicsData::getInstance().getData <std::vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <std::vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <std::vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <std::vector<double> > ("GAMS");
  typemol = 
    PhysicsData::getInstance().getData <std::vector<std::string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
