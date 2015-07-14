// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/ChangeBndryPtr.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Binsrc.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ChangeBndryPtr, Remeshing>
changeBndryPtrProv("ChangeBndryPtr");

//--------------------------------------------------------------------------//

ChangeBndryPtr::ChangeBndryPtr(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

ChangeBndryPtr::~ChangeBndryPtr()
{
}

//--------------------------------------------------------------------------//

void ChangeBndryPtr::setup()
{
  LogToScreen (VERBOSE, "ChangeBndryPtr::setup => start\n");

  LogToScreen (VERBOSE, "ChangeBndryPtr::setup => end\n");
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::unsetup()
{
  LogToScreen (VERBOSE, "ChangeBndryPtr::unsetup()\n");
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::remesh()
{
  LogToScreen (INFO, "ChangeBndryPtr::remesh() \n");

  logfile.Open(getClassName());

  setMeshData();

  setAddress();

  // clear the boolean inside the vector to reset its elements
  // to "false"
  alreadyRemoved.clear();
  alreadyRemoved.resize(nbfac->at(0),false);

  // resize the vector and set its elements to zero
  inactiveFaceInfo.resize(3,nbfac->at(0)+1);

  logfile("Subr ChangeBndryPtr; NBFAC was = ",nbfac->at(0),"\n");
  logfile("Subr ChangeBndryPtr; NBPOIN was = ",nbpoin->at(0),"\n");
  logfile("Subr ChangeBndryPtr; NPOIN was = ",npoin->at(0),"\n");

  for (int IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {

   if (nodcod->at(IPOIN)==-2) {

    // search the i-node IPOIN in inodptr vector if the node has been
    // dis-actived (nodcod=-2)
    lookForNode(IPOIN);

    // remove the i-node
    unsigned ipoin = IPOIN+1; //c++indeces start from 0
    logfile("Removing node", ipoin, "belongs to edges: ");
    logfile((*nodptr)(ipos,1) , (*nodptr)(ipos,2), "\n");
    removeNode(IPOIN);

    // update the i-face=nodptr(ipos,2)
    updateIface(IPOIN);

    // remove the i-face=nodptr(ipos,3)
    removeIface();
   }
  }

  logfile("Subr ChangeBndryPtr; NBFAC is now = ",nbfac->at(0),"\n\n");

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::lookForNode(int IPOIN)
{
  vector <int> iwork_nodptr;
  iwork_nodptr.resize(nbpoin->at(0));

  for (unsigned i=0; i < nbpoin->at(0); i++) {
   iwork_nodptr.at(i) = (*nodptr)(i,0);
  }

  unsigned ipoin = IPOIN+1; //c++ indeces start from 0

  Binsrc findIpoin(ipoin,iwork_nodptr);
  ipos = findIpoin.callBinsrc();
  if (ipos==-1) {
   cout << "ChangeBndryPtr::error => look at ChangeBndryPtr.log file\n";
   logfile("Entry NOT found for ", ipoin, "\n");
   logfile("Entry NOT found for ", ipoin, "\n");
   exit(1);
  }
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::removeNode(int IPOIN)
{
  inode.resize(2);
  for (unsigned K=0; K<2; K++) {
   iface = (*nodptr)(ipos,K+1);

   if ((*bndfac)(0,iface-1) != IPOIN+1) // c++ indeces start from 0
      { inode.at(K) = (*bndfac)(0,iface-1); }
   else { inode.at(K) = (*bndfac)(1,iface-1); } // c++ indeces start from 0
  }
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::updateIface(int IPOIN)
{
  iface = (*nodptr)(ipos,1);

  logfile("Face ",iface, ": ");
  logfile((*bndfac)(0,iface-1)," ",(*bndfac)(1,iface-1),"\n");

  // update one of the two edges of the current face
  (*bndfac)(0,iface-1)=inode.at(0); 

  // if the closer face has been already removed (it happens when there are two 
  // neighbouring phantom nodes) use the edge point of face next to the removed one
  if(alreadyRemoved.at((*nodptr)(ipos,2))) {

   if(inactiveFaceInfo(0,(*nodptr)(ipos,2))!=IPOIN+1) {
    (*bndfac)(1,iface-1)=inactiveFaceInfo(0,(*nodptr)(ipos,2));
   }
   else {(*bndfac)(1,iface-1)=inactiveFaceInfo(1,(*nodptr)(ipos,2));}

   // make in-active the face used to update the edge point
   (*bndfac)(2,inactiveFaceInfo(2,(*nodptr)(ipos,2))-1) = -(*bndfac)(2,inactiveFaceInfo(2,(*nodptr)(ipos,2))-1);
  }

  // else use the edge point of the closer face
//  else { (*bndfac)(0,iface-1)=inode.at(1); }
  else { (*bndfac)(1,iface-1)=inode.at(1); }

  logfile("has been updated with ",(*bndfac)(0,iface-1)," ",(*bndfac)(1,iface-1),"\n");
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::removeIface()
{
  iface = (*nodptr)(ipos,2);

  (*bndfac)(0,iface-1) = (*bndfac)(0,iface-1); // c++ indeces start from 0
  (*bndfac)(1,iface-1) = (*bndfac)(1,iface-1); // c++ indeces start from 0

  if(!alreadyRemoved.at(iface)) {

   logfile("Face ", iface, ": ");
   logfile((*bndfac)(0,iface-1)," ",(*bndfac)(1,iface-1),"\n");
   logfile("has been removed\n");

   (*bndfac)(2,iface-1) = -(*bndfac)(2,iface-1); // c++ indeces start from 0
   alreadyRemoved.at(iface)= true;
   // store infos on the face closer to the removed one
   // it is useful in case of neighbouring phantom nodes
   // @param inactiveFaceInfo(0,iface)   one of the two edge points
   // @param inactiveFaceInfo(1,iface)   the second edge point
   // @param inactiveFaceInfo(2,iface)   ID face
   inactiveFaceInfo(0,iface)=(*bndfac)(0,(*nodptr)(ipos,1)-1);
   inactiveFaceInfo(1,iface)=(*bndfac)(1,(*nodptr)(ipos,1)-1);
   inactiveFaceInfo(2,iface)=(*nodptr)(ipos,1);
  }

  else { logfile("Warning (!) Face: ", iface, ": ");
         logfile((*bndfac)(0,iface-1)," ",(*bndfac)(1,iface-1),"\n");
         logfile("has been already removed\n"); }
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::setAddress()
{
  unsigned start = 0;
  unsigned totsize = nbfac->at(0) + 2 *
                                    PhysicsInfo::getnbShMax() *
                                    PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  nodptr = new Array2D<int> (nbpoin->at(0),3, &nodptrVect->at(start));
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::freeArray()
{
  delete nodptr; delete bndfac;
}

//----------------------------------------------------------------------------//

void ChangeBndryPtr::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbpoin = MeshData::getInstance().getData <vector<unsigned> > ("NBPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  nodptrVect = MeshData::getInstance().getData <vector<int> >("NODPTR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting



