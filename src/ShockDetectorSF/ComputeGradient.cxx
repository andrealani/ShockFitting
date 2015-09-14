// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/ComputeGradient.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

ComputeGradient::ComputeGradient()
{
}

//--------------------------------------------------------------------------//

ComputeGradient::~ComputeGradient()
{

}

//--------------------------------------------------------------------------//

Array2D<double> ComputeGradient::getGrad(vector<double>* primData)
{
  cout << "[P0] [INFO] ComputeGradient::compute()\n";

  setMeshData();

  // assign starting pointers to array
  setAddress();

  Array2D<double> grad(2,npoin->at(0));
  vector<double> gradMod(npoin->at(0));

  // define file plotting gradient distribution
  FILE* plotGradient;
  plotGradient = fopen("gradientDistribution.plt","w");

  fprintf(plotGradient,"%s","TITLE      = Gradient Distribution\n");
  fprintf(plotGradient,"%s","VARIABLES  =  \"x0\" \"x1\" \"GradX\" \"GradY\" \"Grad\"\n");
  fprintf(plotGradient,"%s %u","ZONE N=",npoin->at(0));
  fprintf(plotGradient,"%s %u"," T= \"Gradient distribution\", E = ", nelem->at(0));
  fprintf(plotGradient,"%s",", F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0\n");

  // the ID number of each median dual cell is the same of the corresponding nodeID
  for(unsigned iMedianCell=0; iMedianCell<npoin->at(0);iMedianCell++) {

   unsigned m_CellNodeID = medianDualCellNode->at(iMedianCell);

   double controlVolume = medianDualCellArea->at(iMedianCell);

   for(unsigned iSurrNode=medianDualCellPtr->at(iMedianCell);
                iSurrNode<medianDualCellPtr->at(iMedianCell+1);
                ++iSurrNode                               ) {

    // the evaluation of the gradient is different between the inner and the boundary
    // points
    if(nodcod->at(m_CellNodeID-1)==0 || 
       (nodcod->at(m_CellNodeID-1)!=0 && 
        nodcod->at(medianDualCellNodes->at(iSurrNode)-1)==0)) {

     grad(0,iMedianCell) += 0.5*(primData->at(m_CellNodeID-1)+
                                 primData->at(medianDualCellNodes->at(iSurrNode)-1))*
                                    ((*XY)(1,medianDualCellNodes->at(iSurrNode)-1)-
                                     (*XY)(1,m_CellNodeID-1));
     grad(1,iMedianCell) += 0.5*(primData->at(m_CellNodeID-1)+
                                    primData->at(medianDualCellNodes->at(iSurrNode)-1))*
                                    ((*XY)(0,m_CellNodeID-1)-
                                     (*XY)(0,medianDualCellNodes->at(iSurrNode)-1));   
    }
    else {

     grad(0,iMedianCell) += (5.*primData->at(m_CellNodeID-1)+
                                primData->at(medianDualCellNodes->at(iSurrNode)-1))/6.0*
                                    0.5*((*XY)(1,medianDualCellNodes->at(iSurrNode)-1)-
                                         (*XY)(1,m_CellNodeID-1));
     grad(1,iMedianCell) += (5.*primData->at(m_CellNodeID-1)+
                               primData->at(medianDualCellNodes->at(iSurrNode)-1))/6.0*
                                    0.5*((*XY)(0,m_CellNodeID-1)-
                                         (*XY)(0,medianDualCellNodes->at(iSurrNode)-1));
    }

    grad(0,iMedianCell) /= controlVolume;
    grad(1,iMedianCell) /= controlVolume;

    gradMod.at(iMedianCell) = sqrt(grad(0,iMedianCell)*grad(0,iMedianCell)+
                                grad(1,iMedianCell)*grad(1,iMedianCell));
   }
  }

  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
   for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
    fprintf(plotGradient,"%26.12E %s",(*XY)(IDIM,IPOIN)," ");
   }
   for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
    fprintf(plotGradient,"%26.12E %s",grad(IDIM,IPOIN)," ");
   }
   fprintf(plotGradient,"%26.12E %s",gradMod.at(IPOIN)," ");
   fprintf(plotGradient,"%s","\n");
  }

  for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
   for(unsigned IVERT=0;IVERT<3;IVERT++) {
   fprintf(plotGradient,"%11u",(*celnod)(IVERT,IELEM));
   }
   fprintf(plotGradient,"%s","\n");
  }


  fclose(plotGradient);

  // de-allocate dynamic array
  freeArray();

  return grad;
}

//--------------------------------------------------------------------------//

Array2D<double> ComputeGradient::getNormalizedGrad(vector<double>* primData)
{
  cout << "[P0] [INFO] ComputeGradient::compute()\n";

  setMeshData();

  // assign starting pointers to array
  setAddress();

  Array2D<double> grad(2,npoin->at(0));

  // define file plotting the gradient distribution
  FILE* plotGradient;
  plotGradient = fopen("gradientDistribution.plt","w");

  fprintf(plotGradient,"%s","TITLE      = Gradient Distribution\n");
  fprintf(plotGradient,"%s","VARIABLES  =  \"x0\" \"x1\" \"GradX\" \"GradY\" \n");
  fprintf(plotGradient,"%s %u","ZONE N=",npoin->at(0));
  fprintf(plotGradient,"%s %u"," T= \"Gradient distribution\", E = ", nelem->at(0));
  fprintf(plotGradient,"%s",", F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0\n");

  for(unsigned iMedianCell=0; iMedianCell<npoin->at(0);iMedianCell++) {

   unsigned m_CellNodeID = medianDualCellNode->at(iMedianCell)-1;
   double controlVolume = medianDualCellArea->at(iMedianCell);

   for(unsigned iSurrNode=medianDualCellPtr->at(iMedianCell);
                iSurrNode<medianDualCellPtr->at(iMedianCell+1);
                ++iSurrNode                               ) {

    // the evaluation of the gradient is different between the inner and the boundary
    // points   
    if(nodcod->at(m_CellNodeID-1)==0) {
     grad(0,iMedianCell) += 0.5*(primData->at(m_CellNodeID-1)+
                                 primData->at(medianDualCellNodes->at(iSurrNode)-1))*
                                    ((*XY)(1,medianDualCellNodes->at(iSurrNode)-1)-
                                     (*XY)(1,m_CellNodeID-1));
     grad(1,iMedianCell) += 0.5*(primData->at(m_CellNodeID-1)+
                                    primData->at(medianDualCellNodes->at(iSurrNode)-1))*
                                    ((*XY)(0,m_CellNodeID-1)-
                                     (*XY)(0,medianDualCellNodes->at(iSurrNode)-1));
    }
    else {
     grad(0,iMedianCell) += (5.*primData->at(m_CellNodeID-1)+
                                primData->at(medianDualCellNodes->at(iSurrNode)-1))/6.0*
                                    0.5*((*XY)(1,medianDualCellNodes->at(iSurrNode)-1)-
                                         (*XY)(1,m_CellNodeID-1));
     grad(1,iMedianCell) += (5.*primData->at(m_CellNodeID-1)+
                               primData->at(medianDualCellNodes->at(iSurrNode)-1))/6.0*
                                    0.5*((*XY)(0,m_CellNodeID-1)-
                                         (*XY)(0,medianDualCellNodes->at(iSurrNode)-1));
    }

    grad(0,iMedianCell) /= controlVolume;
    grad(1,iMedianCell) /= controlVolume;

    double gradNorm = sqrt(grad(0,iMedianCell)*grad(0,iMedianCell)+
                           grad(1,iMedianCell)*grad(1,iMedianCell));

    grad(0,iMedianCell) /= gradNorm;
    grad(1,iMedianCell) /= gradNorm;

   }

   for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
    fprintf(plotGradient,"%26.12E %s",(*XY)(IDIM,m_CellNodeID)," ");
   }
   for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
    fprintf(plotGradient,"%26.12E %s",grad(IDIM,iMedianCell)," ");
   }
   fprintf(plotGradient,"%s","\n");
  }

  fclose(plotGradient);

  // de-allocate dynamic array
  freeArray();

  return grad;
}

//--------------------------------------------------------------------------//

void ComputeGradient::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                        PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double>(PhysicsInfo::getnbDim() ,
                           totsize, &coorVect->at(0) );
  celnod = new Array2D<int>((*nvt),
                            nelem->at(0), &celnodVect->at(0) );
}

//--------------------------------------------------------------------------//

void ComputeGradient::freeArray()
{
  delete XY; delete celnod;
}

//--------------------------------------------------------------------------//

void ComputeGradient::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> >("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> >("NELEM");
  medianDualCellArea =
   MeshData::getInstance().getData <vector<double> >("MedianDualCellArea");
  medianDualCellNode =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNode");
  medianDualCellNodes =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNodes");
  medianDualCellPtr =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellPtr");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
