// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/BndryNodePtrFreez.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
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
ObjectProvider<BndryNodePtrFreez, Remeshing> readBndryNodePtrFreezProv("BndryNodePtrFreez");

//--------------------------------------------------------------------------//

BndryNodePtrFreez::BndryNodePtrFreez(const std::string& objectName) :
 Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

BndryNodePtrFreez::~BndryNodePtrFreez()
{
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::setup()
{
  LogToScreen(VERBOSE, "BndryNodePtrFreez::setup() => start\n");

  LogToScreen(VERBOSE, "BndryNodePtrFreez::setup() => end\n");
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::unsetup()
{
  LogToScreen(VERBOSE, "BndryNodePtrFreez::unsetup()\n");
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::remesh()
{
  LogToScreen(INFO, "BndryNodePtrFreez::remesh()\n");

  logfile.Open(getClassName());

  setMeshData();

  setBndryNodePtrFreez();

  // de.allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::setBndryNodePtrFreez()
{
  getnbFreezedPoints();
  getnbBndryPoints();

  nodptrVect->resize(nbpoin->at(0) * 3);
  iwork_.resize(nbpoin->at(0) * 2);
  iworkRank_.resize(nbpoin->at(0));

  // assign arrays used in BndryNodePtrFreez to MeshData 
  setAddress();

  // return vector nodptr
  myroutine();
}                  

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::myroutine()
{
  unsigned IPOIN, last, ifail;
  int ipos;
  vector <int> iwork_nodptr;
  iwork_nodptr.resize(nbpoin->at(0));

  last=0;
  for (IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if ((*nodcod)[IPOIN] > 0 && (*nodcod)[IPOIN]!=999) {
    iwork_.at(last) = IPOIN+1; // c++ vector index starts from 0
    last++;
   }
  }
  logfile("Found ",last," boundary points\n");

  if (last==nbpoin->at(0)) {
  logfile("Found ",nbpoin->at(0)," boundary points\n");
  }
  else {
   cout << "BndryNodePtrFreez::error => LAST =! NBPOIN\n";
   exit(1);
  }

  // iwork_(1:NBPOIN) stores the NBPOIN node numbers
  Isortrx I(iwork_,&nbpoin->at(0));
  iworkRank_ = I.callIsortrx();


  for (IPOIN=0; IPOIN < nbpoin->at(0); IPOIN++) {
   (*nodptr)(IPOIN,0) = iwork_.at(iworkRank_.at(IPOIN));
   iwork_nodptr.at(IPOIN) = (*nodptr)(IPOIN,0);
  }

  ifail = 0;
  unsigned IFACE=0;
  while (IFACE< nbfac->at(0)) {
  unsigned j=0;
  five:
     while (j<2) {
      IPOIN = (*bndfac)(j, IFACE);
      Binsrc B(IPOIN, iwork_nodptr);
      ipos = B.callBinsrc();
      if (ipos ==-1) {
       cout << "BndryNodePtrFreez::error => entry NOT found for " << IPOIN << endl;
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
     if (ifail !=0) { 
      cout << "BndryNodePtrFreez::error => something wrong in CheckBndryPntr\n";
      exit(1);
     }

  return;
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::getnbFreezedPoints()
{
  *nfpoin=0;
  for (unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if((*nodcod)[IPOIN]==999) { (*nfpoin)++;}
  }

  logfile("Found ",(*nfpoin)," freezed points\n");
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::getnbBndryPoints()
{
  nbpoin->at(0)=0;
  for (unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if((*nodcod)[IPOIN]>0) {nbpoin->at(0)++;}
  }
  nbpoin->at(0) = nbpoin->at(0) - (*nfpoin);

  logfile("Found ",nbpoin->at(0)," boundary points\n");
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::freeArray()
{
  delete bndfac; delete nodptr;
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::setAddress()
{
  unsigned start = 0;
  unsigned totsize = nbfac->at(0) + 2 *
                                    PhysicsInfo::getnbShMax() *
                                    PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (3,totsize,&bndfacVect->at(start));
  nodptr = new Array2D<int> (nbpoin->at(0),3, &nodptrVect->at(start));
}

//--------------------------------------------------------------------------//

void BndryNodePtrFreez::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbpoin = MeshData::getInstance().getData <vector<unsigned> > ("NBPOIN");
  nfpoin = MeshData::getInstance().getData <unsigned> ("NFPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  nodptrVect = MeshData::getInstance().getData <vector<int> >("NODPTR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting


