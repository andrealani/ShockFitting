// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/CoNorm.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

CoNorm::CoNorm(const std::string& objectName) : 
  Remeshing(objectName)
{
}

//----------------------------------------------------------------------------//

CoNorm::~CoNorm()
{
  delete ZRoeShd;
}

//----------------------------------------------------------------------------//

void CoNorm::configure(OptionMap& cmap, const std::string& prefix)
{
  Remeshing::configure(cmap, prefix);
}

//----------------------------------------------------------------------------//

void CoNorm::setAddress()
{
  unsigned start = npoin->at(0) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  ZRoeShd =
    new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroe->at(start));
}

//----------------------------------------------------------------------------//

void CoNorm::setSize()
{
  vShNor->resize((*ndim),(*npshmax),(*nshmax));
}

//----------------------------------------------------------------------------//

void CoNorm::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroe = MeshData::getInstance().getData <vector<double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void CoNorm::setPhysicsData()
{
  gam = PhysicsData::getInstance().getData <double> ("GAM");
  gm1 = PhysicsData::getInstance().getData <double> ("GM1");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  model = PhysicsData::getInstance().getData <vector<string> > ("MODEL");
  mixture = PhysicsData::getInstance().getData <vector<string> > ("MIXTURE");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  Rs = PhysicsData::getInstance().getData <vector<double> > ("RS");
  ie = PhysicsData::getInstance().getData <unsigned> ("IE");
  iev = PhysicsData::getInstance().getData <unsigned> ("IEV");
  ix = PhysicsData::getInstance().getData <unsigned> ("IX");
  iy = PhysicsData::getInstance().getData <unsigned> ("IY");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
      PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nSpecPoints =
      PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  typeSh =
      PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  typeSpecPoints =
      PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  XYSh = PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
  vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  SHinSPPs =
      PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

