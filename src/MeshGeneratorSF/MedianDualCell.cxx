// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/MedianDualCell.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MedianDualCell, MeshGenerator>
medianDualCellProv("MedianDualCell");

//--------------------------------------------------------------------------//

MedianDualCell::MedianDualCell(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

MedianDualCell::~MedianDualCell()
{
}

//--------------------------------------------------------------------------//

void MedianDualCell::setup()
{
  LogToScreen(VERBOSE,"MedianDualCell::setup() => start\n");

  LogToScreen(VERBOSE,"MedianDualCell::setup() => end\n");
}

//--------------------------------------------------------------------------//

void MedianDualCell::unsetup()
{
  LogToScreen(VERBOSE,"MedianDualCell::unsetup()\n");
}

//--------------------------------------------------------------------------// 

void MedianDualCell::generate()
{
  LogToScreen(INFO,"MedianDualCell::generate()\n");

  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // open the log file
  logfile.Open(getClassName());

  unsigned ID_1=0; unsigned ID_2=0; unsigned ID_3=0;
  unsigned m_ID=0; unsigned first_surrNodeID=0; unsigned second_surrNodeID=0;

  // working vector storing element centroid coordinates
  vector<double> centroid(2);
  // working vector storing element edge middle points coordinates
  vector<double> middleEdge_1(2);
  vector<double> middleEdge_2(2);

  // the boundary points are not taken into account
  medianDualCellPtr->resize(npoin->at(0)+1);
  medianDualCellArea->resize(npoin->at(0));
  medianDualCellNode->resize(npoin->at(0));

  unsigned inPoin = 0;
  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {

   unsigned m_poin = IPOIN+1;
   unsigned ptrCounter = 0;

   // working vector storing nodes of a median dual cell
   vector<unsigned> surrNodes; //(1,0);

    logfile("Building median dual cell for nodeID:: ",m_poin,"\n");

    for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {

     for(unsigned IVERT=0;IVERT<(*nvt);IVERT++) {

      if(m_poin==(*celnod)(IVERT,IELEM)) {

       if(IVERT==0) {ID_1=0; ID_2 = 1; ID_3=2;
                     m_ID = IVERT; first_surrNodeID = ID_2; second_surrNodeID=ID_3;}
       if(IVERT==1) {ID_1=IVERT-1; ID_2=IVERT; ID_3=IVERT+1;
                     m_ID = IVERT; first_surrNodeID = ID_1; second_surrNodeID=ID_3;}
       if(IVERT==2) {ID_1=IVERT-2; ID_2=IVERT-1; ID_3=IVERT;
                     m_ID = IVERT; first_surrNodeID = ID_1; second_surrNodeID=ID_2;}
       logfile("..........................");
       logfile("Cell: ",IELEM+1, " => nodes: ");
       logfile((*celnod)(ID_1,IELEM)," ",
               (*celnod)(ID_2,IELEM)," ",
               (*celnod)(ID_3,IELEM),"\n"); 

       centroid.at(0) = ((*XY)(0,(*celnod)(ID_1,IELEM)-1)+
                         (*XY)(0,(*celnod)(ID_2,IELEM)-1)+
                         (*XY)(0,(*celnod)(ID_3,IELEM)-1))*1/3;
       centroid.at(1) = ((*XY)(1,(*celnod)(ID_1,IELEM)-1)+
                         (*XY)(1,(*celnod)(ID_2,IELEM)-1)+
                         (*XY)(1,(*celnod)(ID_3,IELEM)-1))*1/3;

       logfile("centroid= (",centroid.at(0),", ",centroid.at(1),")\n");
        
       middleEdge_1.at(0) = 0.5*((*XY)(0,(*celnod)(IVERT,IELEM)-1)+
                                 (*XY)(0,(*celnod)(first_surrNodeID,IELEM)-1));
       middleEdge_1.at(1) = 0.5*((*XY)(1,(*celnod)(IVERT,IELEM)-1)+
                                 (*XY)(1,(*celnod)(first_surrNodeID,IELEM)-1));
       middleEdge_2.at(0) = 0.5*((*XY)(0,(*celnod)(IVERT,IELEM)-1)+
                                 (*XY)(0,(*celnod)(second_surrNodeID,IELEM)-1));
       middleEdge_2.at(1) = 0.5*((*XY)(1,(*celnod)(IVERT,IELEM)-1)+
                                (*XY)(1,(*celnod)(second_surrNodeID,IELEM)-1));
      
       double myarea = 0.5*(((*XY)(0,(*celnod)(m_ID,IELEM)-1)-centroid.at(0))*
                             (middleEdge_2.at(1)-middleEdge_1.at(1))+
                            (middleEdge_1.at(0)-middleEdge_2.at(0))*
                             ((*XY)(1,(*celnod)(m_ID,IELEM)-1)-centroid.at(1)));
       if (myarea<0) { 
        myarea = 0.5*(((*XY)(0,(*celnod)(m_ID,IELEM)-1)-centroid.at(0))*
                       (middleEdge_1.at(1)-middleEdge_2.at(1))+
                      (middleEdge_2.at(0)-middleEdge_1.at(0))*
                       ((*XY)(1,(*celnod)(m_ID,IELEM)-1)-centroid.at(1)));
       }

       medianDualCellArea->at(inPoin) += myarea;

       if(medianDualCellArea->at(inPoin)<0) {
        cout << "MedianDualCell::error => median dual cell area of nodeID ";
        cout << m_poin << " is negative\n";
        cout << "                         check if the element nodes are \n";
        cout << "                         numbered in a counterclockwise order\n";
        cout << "                         inside the MedialDualCell log file\n";
        exit(1);
       }

       surrNodes.insert(surrNodes.begin()+ptrCounter,
                        (*celnod)(first_surrNodeID,IELEM));
       surrNodes.insert(surrNodes.begin()+ptrCounter+1,
                        (*celnod)(second_surrNodeID,IELEM));

       ptrCounter+=2;

      } // if m_poin==celnod   
     } //for IVERT=0:nvt
    } // for IELEM=0:nelem

    logfile("\nControl volume of the median dual cell= ",
            medianDualCellArea->at(inPoin),"\n");

    medianDualCellNode->at(inPoin) = m_poin;

    eraseDuplicate(surrNodes);

    medianDualCellPtr->at(inPoin+1) = medianDualCellPtr->at(inPoin)+
                                      surrNodes.size();

    medianDualCellNodes->resize(medianDualCellPtr->at(inPoin)+surrNodes.size());
 
    logfile("The nodes of the median dual cell are => ");
    // store the nodes building the median dual cell
    for(unsigned iSurrNode=0;iSurrNode<surrNodes.size();iSurrNode++) {
     medianDualCellNodes->at(medianDualCellPtr->at(inPoin)+
                            iSurrNode)=surrNodes.at(iSurrNode);
     logfile(surrNodes.at(iSurrNode)," ");
    }
 
    logfile("\n--------------------------------------------\n\n");

    inPoin++;

  } // for IPOIN=0:npoin

  // resize the medianDualCellNodes vector 
  medianDualCellNodes->resize(medianDualCellPtr->at(npoin->at(0)));

  logfile("Number of median dual cells is ",medianDualCellNode->size());

  // close the log file
  logfile.Close();

  if(medianDualCellNode->size()!=npoin->at(0)) {
   cout << "MedianDualCell::error => number of evaluated median dual cells\n";
   cout << "                         is not equal to the number of points\n";
   cout << "                         check the MedianDualCell log file\n";
   exit(1);
  }

  // de-allocate the dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void MedianDualCell::eraseDuplicate(vector<unsigned>& m_vector)
{
  vector<unsigned> m_vectorNoDuplicate(m_vector.size());

  unsigned counter=0; unsigned duplicate=0;
  for(unsigned j=0;j<m_vector.size();j++){
   for(unsigned k=j+1;k<m_vector.size();k++) {
    if(m_vector.at(j)==m_vector.at(k)) { duplicate++; }
   }
   if(duplicate==0) { m_vectorNoDuplicate.at(counter)=m_vector.at(j);
                      counter++;
                    }
   else { duplicate=0; }
  }

  m_vector.resize(counter);
  m_vectorNoDuplicate.resize(counter);

  for(unsigned i=0;i<m_vectorNoDuplicate.size();i++) {
   m_vector.at(i) = m_vectorNoDuplicate.at(i);
  }
}

//--------------------------------------------------------------------------//

void MedianDualCell::generate(std::string)
{
}

//--------------------------------------------------------------------------//

void MedianDualCell::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                        PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double>(PhysicsInfo::getnbDim() ,
                           totsize, &coorVect->at(0) ); 
  celnod = new Array2D<int> ((*nvt), nelem->at(0), &celnodVect->at(0));
}

//--------------------------------------------------------------------------//

void MedianDualCell::freeArray()
{
  delete XY; delete celnod;
}

//--------------------------------------------------------------------------//

void MedianDualCell::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  medianDualCellArea = 
   MeshData::getInstance().getData <vector<double> >("MedianDualCellArea");
  medianDualCellNode =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNode");
  medianDualCellNodes =
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellNodes");
  medianDualCellPtr = 
   MeshData::getInstance().getData <vector<unsigned> >("MedianDualCellPtr");
}

//--------------------------------------------------------------------------//

void MedianDualCell::setPhysicsData()
{
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
