// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/BndryNodePtr.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Isortrx.hh"
#include "MathTools/Binsrc.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<BndryNodePtr, Remeshing> readBndryNodePtrProv("BndryNodePtr");

//--------------------------------------------------------------------------//

BndryNodePtr::BndryNodePtr(const std::string& objectName)
             :Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

BndryNodePtr::~BndryNodePtr()
{
}

//--------------------------------------------------------------------------//

void BndryNodePtr::setup()
{
  LogToScreen(VERBOSE, "BndryNodePtr::setup() => start\n");

  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "BndryNodePtr::setup() => end\n");
}

//--------------------------------------------------------------------------//

void BndryNodePtr::unsetup()
{
  LogToScreen(VERBOSE, "BndryNodePtr::unsetup()\n");

  logfile.Close();
}

//--------------------------------------------------------------------------//

void BndryNodePtr::remesh()
{
  LogToScreen(INFO, "BndryNodePtr::remesh()\n");

  setMeshData();

  setBndryNodePtr();
}

//--------------------------------------------------------------------------//

void BndryNodePtr::setBndryNodePtr()
{
  getnbFreezedPoints();
  getnbBndryPoints();
  
  nodptr->resize((*nbpoin),3);
  iwork_.resize((*nbpoin));
  iworkRank_.resize((*nbpoin));

  // return vector nodptr
  myroutine();
}                  

//--------------------------------------------------------------------------//

void BndryNodePtr::myroutine()
{
  unsigned IPOIN, last, ifail;
  int ipos;
  vector <int> iwork_nodptr;
  iwork_nodptr.resize((*nbpoin));

  last=0;
  for (IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   if ((*nodcod)[IPOIN] > 0 && (*nodcod)[IPOIN]!=999) {
    iwork_.at(last) = IPOIN+1; // c++ vector index starts from 0
    last++;
   }
  }
  logfile("Found ",last," boundary points\n");

  if (last==(*nbpoin)) {
  logfile("Found ",(*nbpoin)," boundary points\n");
  }
  else {
   cout << "LAST =! NBPOIN\n";
   exit(1);
  }

  // iwork_(1:NBPOIN) stores the NBPOIN node numbers
  Isortrx I(iwork_,nbpoin);
  iworkRank_ = I.callIsortrx();

  for (IPOIN=0; IPOIN < (*nbpoin); IPOIN++) {
   (*nodptr)(IPOIN,0) = iwork_.at(iworkRank_.at(IPOIN));
   iwork_nodptr.at(IPOIN) = (*nodptr)(IPOIN,0);
  }

  ifail = 0;
  unsigned IFACE=0;
  while (IFACE< (*nbfac)) {
  unsigned j=0;
  five:
     while (j<2) {
      IPOIN = (*bndfac)(j, IFACE);
      Binsrc B(IPOIN, iwork_nodptr);
      ipos = B.callBinsrc();
      if (ipos ==-1) {
       cout << "Subr. CheckBndryPntr: Entry NOT found for " << IPOIN << endl;
       cout << "Color: " << (*bndfac)(2,IFACE) << endl;
       cout << (*bndfac)(0,IFACE) << " , " << (*bndfac)(1,IFACE) << endl;
       exit(1);
      }
      for (unsigned k=1; k<3; k++) {
       if ((*nodptr)(ipos,k)==0) {
        (*nodptr)(ipos,k)=IFACE+1;// c++ indeces start from 0
        j++;
        goto five;
        break;
       }
      }
  // node seems to belong to two boundary faces
      ifail = IPOIN;
      cout << "Node seems to belong to two boundary faces\n";
      cout << "Face no. " << IFACE << endl;
      for (unsigned k=0; k<3; k++) { cout << (*bndfac)(k,IFACE) << ", " << endl;}
      cout << "Node no. " << IPOIN << endl;
      for (unsigned k=0; k<3; k++) { cout << (*nodptr)(ipos,k) << ", " << endl;}
      cout << "Color " << (*bndfac)(2,IFACE) << endl;
      goto seven;
     }
     IFACE++;
  }

  seven:
     if (ifail !=0) { cout << "Unrecoverable error in CheckBndryPntr\n"; exit(1);}

  return;
}

//--------------------------------------------------------------------------//

void BndryNodePtr::getnbFreezedPoints()
{
  *nfpoin=0;
  for (unsigned IPOIN=0; IPOIN<(*npoin); IPOIN++) {
   if((*nodcod)[IPOIN]==999) { (*nfpoin)++;}
  }
}

//--------------------------------------------------------------------------//

void BndryNodePtr::getnbBndryPoints()
{
  *nbpoin=0;
  for (unsigned IPOIN=0; IPOIN<(*npoin)-1; IPOIN++) {
   if((*nodcod)[IPOIN]>0) {(*nbpoin)++;}
  }
  (*nbpoin) = (*nbpoin) - (*nfpoin);
}

//--------------------------------------------------------------------------//

void BndryNodePtr::setMeshData()
{
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  nbpoin = MeshData::getInstance().getData <unsigned> ("NBPOIN");
  nfpoin = MeshData::getInstance().getData <unsigned> ("NFPOIN");
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  nodcod = MeshData::getInstance().getData< std::vector<int> >("NODCOD");
  nodptr = MeshData::getInstance().getData< Array2D<int> >("NODPTR");
  bndfac = MeshData::getInstance().getData< Array2D<int> >("BNDFAC");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting



