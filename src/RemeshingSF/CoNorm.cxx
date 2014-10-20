// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/CoNorm.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
//ObjectProvider<CoNorm, Remeshing> computeNormProv("CoNorm");

//--------------------------------------------------------------------------//

CoNorm::CoNorm(const std::string& objectName)
 :NormalUnitVect(objectName)
{

}

//----------------------------------------------------------------------------//

CoNorm::~CoNorm()
{
}

//----------------------------------------------------------------------------//

void CoNorm::setAddress()
{
  unsigned start = (*npoin) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  r_ZRoeShd =
    new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroe->at(start));
}

//----------------------------------------------------------------------------//

void CoNorm::setSize()
{
  r_vShNor->resize((*ndim),(*npshmax),(*nshmax));
}

//----------------------------------------------------------------------------//

void CoNorm::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
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
  r_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  r_nShockPoints =
      PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  r_nSpecPoints =
      PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  r_typeSh =
      PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  r_typeSpecPoints =
      PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  r_XYSh = PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
  r_vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  r_SHinSPPs =
      PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

