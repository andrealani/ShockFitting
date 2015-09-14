// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/DensityGradientMaximaPg2D.hh"
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

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DensityGradientMaximaPg2D, ShockDetector>
densityGradientMaxPg2dProv("DensityGradientMaximaPg2D");

//--------------------------------------------------------------------------//

DensityGradientMaximaPg2D::DensityGradientMaximaPg2D(const std::string& objectName) :
  ShockDetector(objectName)
{
  m_filteringFactor = 0;
  addOption("filteringFactor",&m_filteringFactor,
             "Value of the chosen filtering factor");
}

//--------------------------------------------------------------------------//

DensityGradientMaximaPg2D::~DensityGradientMaximaPg2D()
{
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::setup()
{
  LogToScreen(VERBOSE,"DensityGradientMaximaPg2D::setup() => start\n");

  LogToScreen(VERBOSE,"DensityGradientMaximaPg2D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::unsetup()
{
  LogToScreen(VERBOSE,"DensityGradientMaximaPg2D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::detect()
{
  LogToScreen(INFO,"DensityGradientMaximaPg2D::detect()\n");
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::detect(std::vector<double>& primitiveVar)
{
  LogToScreen(INFO,"DensityGradientMaximaPg2D::detect()\n");

  setMeshData();
  setPhysicsData();

  // assign starting pointers to the array
  setAddress();

  // resize vectors and array
  setSize();

  setPressureID();
  setVelocityID();
  setTemperatureID();

  // create object evaluating density gradient
  ComputeGradient gradRHO;
  ComputeGradient gradFirstRHO;

  // create the array storing the density gradients in the median dual cells
  Array2D<double> densityGradient;
  Array2D<double> firstDensityDerivativeGradient;

  vector<double> firstDensityDerivative(npoin->at(0));
  vector<double> secondDensityDerivative(npoin->at(0));

  // create working vectors storing pimitive components
  vector<double> density(npoin->at(0));
  vector<double> u(npoin->at(0));
  vector<double> v(npoin->at(0));

  /// the algorithm stores the shock points distribution as if only one
  /// shock is detected. The number of shock will be evaluated after 
  (*nShocks) = 1;
  unsigned ISH = 0;

  unsigned iShPoin=0;

  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
   density.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+pID)/
                       (ReferenceInfo::getRgas()*
                        primitiveVar.at(IPOIN*(*ndof)+tID));
   u.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+uID);
   v.at(IPOIN) = primitiveVar.at(IPOIN*(*ndof)+vID);
  }

  // compute density gradients
  densityGradient = gradRHO.getGrad(&density);

  for(unsigned iMedianCell=0; iMedianCell<npoin->at(0);iMedianCell++) {

   unsigned m_CellNodeID = medianDualCellNode->at(iMedianCell)-1;

   double velocityMod = sqrt(u.at(m_CellNodeID)*u.at(m_CellNodeID)+
                             v.at(m_CellNodeID)*v.at(m_CellNodeID));
   u.at(m_CellNodeID) /= velocityMod;
   v.at(m_CellNodeID) /= velocityMod;
   firstDensityDerivative.at(iMedianCell) = densityGradient(0,iMedianCell)*
                                            u.at(m_CellNodeID)+
                                            densityGradient(1,iMedianCell)*
                                            v.at(m_CellNodeID);
  }

  firstDensityDerivativeGradient = gradFirstRHO.getGrad(&firstDensityDerivative);

  for(unsigned iMedianCell=0; iMedianCell<npoin->at(0);iMedianCell++) {

   unsigned m_CellNodeID = medianDualCellNode->at(iMedianCell)-1;

   secondDensityDerivative.at(iMedianCell) = 
    firstDensityDerivativeGradient(0,iMedianCell)*u.at(m_CellNodeID)+
    firstDensityDerivativeGradient(1,iMedianCell)*v.at(m_CellNodeID);

   // check if the condition on the density derivatives is satisfied
   if(abs(secondDensityDerivative.at(iMedianCell))<1.e-3 &&
      firstDensityDerivative.at(iMedianCell)>m_filteringFactor) {
    (*XYSh)(0,iShPoin,ISH) = (*XY)(0,m_CellNodeID);
    (*XYSh)(1,iShPoin,ISH) = (*XY)(1,m_CellNodeID);
    iShPoin++;
   }
  }

  nShockPoints->at(ISH) = iShPoin;
  nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::setSize()
{
  nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  XYSh->resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0)); 
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::freeArray()
{
  delete XY;
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> >("NPOIN");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  medianDualCellNode =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNode");
  medianDualCellPtr =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellPtr");
}

//--------------------------------------------------------------------------//

void DensityGradientMaximaPg2D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned>("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
