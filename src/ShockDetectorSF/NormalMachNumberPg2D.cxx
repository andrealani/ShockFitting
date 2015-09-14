// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/NormalMachNumberPg2D.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "MathTools/LeastSquaresMethod.hh"
#include "SConfig/ObjectProvider.hh"
#include "ShockDetectorSF/ComputeGradient.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------/

// this variable instantiation activates the self-registration mechanism
ObjectProvider<NormalMachNumberPg2D, ShockDetector>
normalMachNumberPgProv("NormalMachNumberPg2D");

//--------------------------------------------------------------------------//

NormalMachNumberPg2D::NormalMachNumberPg2D(const std::string& objectName) :
  ShockDetector(objectName)
{
}

//--------------------------------------------------------------------------/

NormalMachNumberPg2D::~NormalMachNumberPg2D()
{
}

//--------------------------------------------------------------------------//

void NormalMachNumberPg2D::setup()
{
  LogToScreen (VERBOSE, "NormalMachNumberPg2D::setup => start\n");

  LogToScreen (VERBOSE, "NormalMachNumberPg2D::setup => end\n");
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::unsetup()
{
  LogToScreen(VERBOSE, "NormalMachNumberPg2D::unsetup()\n");
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::detect()
{
  LogToScreen(INFO,"NormalMachNumberPg2D::detect()\n");
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::detect(vector<double>& primitiveVar)
{
  LogToScreen(INFO,"NormalMachNumberPg2D::detect()\n");

  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // resize vectors and array
  setSize();

  setPressureID();
  setVelocityID();
  setTemperatureID();

  // create the object evaluating pressure gradient
  ComputeGradient gradP;

  // create the array storing the pressure gradients in the median dual cells
  Array2D<double> pressureGradient;

  // create working vectors storing velocity components and pressure
  vector<double> pressure(npoin->at(0));
  vector<double> u(npoin->at(0));
  vector<double> v(npoin->at(0));
  vector<double> a(npoin->at(0));

  vector<double> normalMach(npoin->at(0));

  /// the algorithm stores the shock points distribution as if only one
  /// shock is detected. The number of shock will be evaluated after 
  (*nShocks) = 1;
  unsigned ISH = 0;

  unsigned iShPoin=0;

  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
   pressure.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+pID);
   u.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+uID);
   v.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+vID);
   a.at(IPOIN) = sqrt(ReferenceInfo::getgam()*
                      ReferenceInfo::getRgas()*
                      primitiveVar.at(IPOIN*(*ndof)+tID));
  }

  // compute pressure gradients
  pressureGradient = gradP.getGrad(&pressure);

  for(unsigned iMedianCell=0; iMedianCell<npoin->at(0);iMedianCell++) {

   unsigned m_CellNodeID = medianDualCellNode->at(iMedianCell)-1;
   double normPressureGrad = sqrt(pow(pressureGradient(0,iMedianCell),2)+
                                  pow(pressureGradient(1,iMedianCell),2));

   /// (!) check if abs must be used
   double Mx = abs(u.at(m_CellNodeID)/a.at(m_CellNodeID));
   double My = abs(v.at(m_CellNodeID)/a.at(m_CellNodeID));

   // compute the normal Mach number from pressure distribution
   normalMach.at(iMedianCell) = Mx*pressureGradient(0,iMedianCell)+
                                My*pressureGradient(1,iMedianCell);
   normalMach.at(iMedianCell) /= normPressureGrad; 

   // check if the shock condition is satisfied
   if(abs(normalMach.at(iMedianCell)-1)<=0.4) {
     (*XYSh)(0,iShPoin,ISH) = (*XY)(0,m_CellNodeID);
     (*XYSh)(1,iShPoin,ISH) = (*XY)(1,m_CellNodeID);
     iShPoin++;
   }
  }

  nShockPoints->at(ISH) = iShPoin;
  nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;

  // de-allocate dynamic array
  freeArray();
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::setSize()
{
  nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  XYSh->resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax());
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0));
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::freeArray()
{
  delete XY;
}

//--------------------------------------------------------------------------//

void NormalMachNumberPg2D::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> >("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> >("NELEM");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  medianDualCellNode =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNode");
  medianDualCellPtr =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellPtr");  
}

//----------------------------------------------------------------------------//

void NormalMachNumberPg2D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned>("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
