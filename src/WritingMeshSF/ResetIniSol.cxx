// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "WritingMeshSF/ResetIniSol.hh"
#include "SConfig/ObjectProvider.hh"
#include "WritingMeshSF/Burgy.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ResetIniSol, WritingMesh> resetIniSolProv("ResetIniSol");

//--------------------------------------------------------------------------//

ResetIniSol::ResetIniSol(const std::string& objectName)
 : WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

ResetIniSol:~ResetIniSol()
{
}

//--------------------------------------------------------------------------//

void ResetIniSol::setup()
{
  LogToScreen(VERBOSE,"ResetIniSol::setup() => start\n");

  LogToScreen(VERBOSE,"ResetIniSol::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ResetIniSol::setup()
{
  LogToScreen(VERBOSE,"ResetIniSol::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ResetIniSol::write()
{
  LogToScreen(INFO,"ResetIniSol::write()\n");

  setPhysicsData();

  setMeshData();

  setAddress();

  unsigned N = 0;

  Burgy B;

  for (unsigned I=0; I<(*npoin); I++) {
   if (std::abs((*w_XY)(1,I)<=10e-8) {
    ++N;
    (*zroe)(0,I) = B.getBurger((*XY)(0,I));
   }
  }
  cout << "Setting exact sol in " << N << " gridpoints" << endl;
  

}

//--------------------------------------------------------------------------//

void ResetIniSol::setAddress()
{
  unsigned start;
  start = (*npoin)*(*ndof);
  r_ZRoeShu = 
    new Array3D <double> ((*ndof),(*npshmax),(*nshmax),&zroe->at(start));
  start = (*npoin) * (*ndof) + (*npshmax) * (*nshmax) * (*ndof);
  r_ZRoeShd = 
    new Array3D <double> ((*ndofmax),(*npshmax),(*nshmax),&zroe->at(start));
  start = 0;
  w_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
  start = (*npoin) * (*ndim);
  w_XYShu =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
  start = (*npoin) * (*ndim) + (*npshmax) * (*nshmax) * (*ndim);
  w_XYShd =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));

}

//--------------------------------------------------------------------------//

void ResetIniSol::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN"); 
  zroe = MeshData::getInstance().getData< std::vector<double> >("ZROE");
  coor = MeshData::getInstance().getData< std::vector<double> >("COOR");
}

//--------------------------------------------------------------------------//

void ResetIniSol::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ndofmax = PhysicsData::getInstance().getData <unsigned> ("NDOFMAX");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  w_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks"); 
  w_nShockPoints = 
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

