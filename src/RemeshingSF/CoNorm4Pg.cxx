// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "RemeshingSF/CoNorm4Pg.hh"
#include "RemeshingSF/ShpDpndnc.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//---------------------------------------------------------------------------//

namespace ShockFitting {

//---------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CoNorm4Pg, CoNorm> computeNormalVector4PgProv("CoNorm4Pg");

//--------------------------------------------------------------------------//

CoNorm4Pg::CoNorm4Pg(const std::string& objectName) :
  CoNorm(objectName)
{
}

//----------------------------------------------------------------------------//

CoNorm4Pg::~CoNorm4Pg()
{
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setup()
{
  LogToScreen(VERBOSE,"CoNorm4Pg::setup() => start\n");

  LogToScreen(VERBOSE,"CoNorm4Pg::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::unsetup()
{
  LogToScreen(VERBOSE,"CoNorm4Pg::unsetup()\n");
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::remesh()
{
  LogToScreen(INFO,"CoNorm4Pg::remesh()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setSize();

  // write status on log file
  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {

    logfile("I: ",(*ZRoeShd)(0,I,ISH), " " ,(*ZRoeShd)(1,I,ISH));
    logfile("I: ",(*ZRoeShd)(2,I,ISH), " " ,(*ZRoeShd)(3,I,ISH));
   }
  }


  // compute normal vector for each shock
  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {

     computeTau(ISH,I);

     // assign normal vector
     (*vShNor)(0,I,ISH) = -tauy;
     (*vShNor)(1,I,ISH) = taux;
   }
  }

  // compute normal vectors for typeSh="S"
  setVShNorForStype();

  // fix normal vectors for special points
  // it forces the direction of the contact discontinuity normal vector
  // in order to define an angle equal to 90° with mach stem
  // normal vector
  for (unsigned ISPPNTS=0; ISPPNTS<(*nSpecPoints); ISPPNTS++) {

   // special point: wall point without reflection
   if (typeSpecPoints->at(ISPPNTS)=="WPNRX") {setVShNorForWPNRX(ISPPNTS);}

   // special point: connection between two shocks
   else if (typeSpecPoints->at(ISPPNTS)=="C") {setVShNorForC(ISPPNTS);}

   // special point: triple point
   else if (typeSpecPoints->at(ISPPNTS)=="TP") {setVShNorForTP(ISPPNTS);}

   // write the computed normal vectors on tecplot file
   writeTecPlotFile();
  }

  logfile.Close();
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::computeTau(unsigned ISH, unsigned I)
{
  unsigned J, J2, ipoin;
  ush = 0; vsh = 0;
  xi = (*XYSh)(0,I,ISH);
  yi = (*XYSh)(1,I,ISH);

  if (I < (nShockPoints->at(ISH)-1)) {
   // one point forward
   J=I+1;
   // coordinates of the one point forward
   onePointForward(J, ISH);

   // recover status for the forward point
   ipoin = J+1; // c++ indeces start from 0
   recoverState("forward",J,ISH);

   if (I < (nShockPoints->at(ISH)-1)) {
    // two points forward
    J2=I+2;
    // coordinates of two points forward
    twoPointsForward(J2,ISH);
   } // J2
   else { setTauIp2ToZero();}
  } // J
  else { // I = nShockPoints(ISH)-1 
   setTauIp1ToZero();
   setTauIp2ToZero();
  }

  // evaluate dependency
  // if I is the last point it is computed in backward direction
  if (I!=0 && I!= nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepip(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depip1 = shockDepip.getDependence();
  }
  else if (I==nShockPoints->at(ISH)) {depip1=0; depim1=1;}

  if (I>0) {
   // one point backward
   J=I-1;
   // coordinates of the one point backward
   onePointBackward(J,ISH);

   // recover status for the backward point
   ipoin = J+1; // c++ indeces start from 0
   recoverState("backward",J,ISH);

   if (I>1) {
    // two points backward
    J2=I-2;
    // coordinates of two points backward
    twoPointsBackward(J2,ISH);
   }
   else {setTauIm2ToZero();}
  }
  else { // I=0
   setTauIm1ToZero();
   setTauIm2ToZero();
  }

  // evaluate dependency
  // if I is the first point it is computed in forward direction
  if (I!=0 && I!= nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepim(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depim1 = shockDepim.getDependence();
  }
  else if (I==0) {depip1=1; depim1=0;}

  setLp();
  setLm();

  if(I==0) {depim1=0; depip1=1; lm12=1;}
  if(I==nShockPoints->at(ISH)-1) {depim1=1; depip1=0; lp12=1.0;}

  taux = (depim1*tauxim1*lp12+depip1*tauxip1*lm12);
  tauy = (depim1*tauyim1*lp12+depip1*tauyip1*lm12);

  ipoin = I+1; // c++ indeces start from 0
  logfile("\n", ipoin, "tau ", taux, "   ", tauy);
  logfile("           ", depim1, " ", depip1, "\n\n");
  logfile ("\n            --------------          \n");
  tau = sqrt(taux*taux+tauy*tauy);
  taux = taux/tau;
  tauy = tauy/tau;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setVShNorForStype()
{
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   if(typeSh->at(ISH)=="S") {
    unsigned ii = 0;
    for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
     ui = (*ZRoeShd)(2,I,ISH)/(*ZRoeShd)(0,I,ISH);
     vi = (*ZRoeShd)(3,I,ISH)/(*ZRoeShd)(0,I,ISH);
     dum = ui*(*vShNor)(0,I,ISH)+vi*(*vShNor)(1,I,ISH);
     if(dum>0) {ii++;}
    }
    if (ii < nShockPoints->at(ISH)/2 - 1) { break; }
     for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
      (*vShNor)(0,I,ISH) = -(*vShNor)(0,I,ISH);
      (*vShNor)(1,I,ISH) = -(*vShNor)(1,I,ISH);
     } // I
   } // if typeSh->at(ISH)=="S"
  } // for
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setVShNorForWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  (*vShNor)(0,IP.at(0),ISH.at(0)) =
      (*vShNor)(0,IP.at(0),ISH.at(0))/abs((*vShNor)(0,IP.at(0),ISH.at(0)));
  (*vShNor)(1,IP.at(0),ISH.at(0)) = 0;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setVShNorForC(unsigned ISPPNTS)
{
  setShockIndeces(2,ISPPNTS);
  
  nx1 = (*vShNor)(0,IP.at(0),ISH.at(0));
  ny1 = (*vShNor)(1,IP.at(0),ISH.at(0));
  nx2 = (*vShNor)(0,IP.at(1),ISH.at(1));
  ny2 = (*vShNor)(1,IP.at(1),ISH.at(1));
  
  nx1 = nx1+nx2;
  ny1 = ny1+ny2;
  
  dum = sqrt(nx1*nx1+ny1*ny1);
  nx1 = nx1/dum;
  ny1 = ny1/dum;
  
  (*vShNor)(0,IP.at(0),ISH.at(0)) = nx1;
  (*vShNor)(1,IP.at(0),ISH.at(0)) = ny1;
  (*vShNor)(0,IP.at(1),ISH.at(0)) = nx1;
  (*vShNor)(1,IP.at(1),ISH.at(0)) = ny1;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setVShNorForTP(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  // ISH.at(2) Mach stem
  // ISH.at(3) contact discontinuity
  setShockIndeces(4,ISPPNTS);

  nx2 = (*vShNor)(0,IP.at(1),ISH.at(1));
  ny2 = (*vShNor)(1,IP.at(1),ISH.at(1));

  nx4 = (*vShNor)(0,IP.at(3),ISH.at(3));
  ny4 = (*vShNor)(1,IP.at(3),ISH.at(3));

  dum = nx2*nx4+ny2*ny4;
  if (dum < 0) {
   for (unsigned I=0; I<nShockPoints->at(ISH.at(3)); I++) {
    (*vShNor)(0,I,ISH.at(3)) = -(*vShNor)(0,I,ISH.at(3));
    (*vShNor)(1,I,ISH.at(3)) = -(*vShNor)(1,I,ISH.at(3));
   }
  }
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::writeTecPlotFile()
{
  ofstream tecfile;
  tecfile.open("shocknor.dat");

  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   tecfile << "TITLE = Shock normals\n";
   tecfile << "VARIABLES = X Y Z(1) Z(2) NX NY\n";
   tecfile << "ZONE T='sampletext', F = FEPOINT, ET = TRIANGLE ";
   tecfile << "N = " << nShockPoints->at(ISH);
   tecfile << ", E = " << nShockPoints->at(ISH)-1 << "\n";
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++)
     {tecfile << (*XYSh)(K,I,ISH) << " ";}
    tecfile << "\n";
    tecfile << 1 << " " << 1 << " ";
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++) 
     {tecfile << (*vShNor)(K,I,ISH) << " ";}
    tecfile << "\n";
   }
   for (unsigned I=0; I<nShockPoints->at(ISH)-1; I++) {
    tecfile << I << " " << I+1 << " " << I << "\n";
   }
  }
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::recoverState(string direction, unsigned J, unsigned ISH)
{
  uj = (*ZRoeShd)(2,J,ISH)/(*ZRoeShd)(0,J,ISH);
  vj = (*ZRoeShd)(3,J,ISH)/(*ZRoeShd)(0,J,ISH);
  roj = (*ZRoeShd)(0,J,ISH)*(*ZRoeShd)(0,J,ISH);
  help = ((*ZRoeShd)(2,J,ISH),2)+pow((*ZRoeShd)(3,J,ISH),2);
  pj = PhysicsInfo::getGm1() / PhysicsInfo::getGam() *
       ((*ZRoeShd)(0,J,ISH)*(*ZRoeShd)(1,J,ISH)-0.5*help);
  aj = sqrt(PhysicsInfo::getGam()*pj/roj);

  if (typeSh->at(ISH)=="D") { aj = 0;}
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
{
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setLp()
{
  lp12 = pow(tauxip1,2)+pow(tauyip1,2);
  lp22 = pow(tauxip2,2)+pow(tauyip2,2);
  lp1 = sqrt(lp12);
  lp2 = sqrt(lp22); 
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setLm()
{
  lm12 = pow(tauxim1,2)+pow(tauyim1,2);
  lm22 = pow(tauxim2,2)+pow(tauyim2,2);
  lm1 = sqrt(lm12);
  lm2 = sqrt(lp22);
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::onePointForward(unsigned J, unsigned ISH)
{
  xj = (*XYSh)(0,J,ISH);
  yj = (*XYSh)(1,J,ISH);
  tauxip1 = xj-xi;
  tauyip1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::twoPointsForward(unsigned J2, unsigned ISH)
{
  xj2 = (*XYSh)(0,J2,ISH);
  yj2 = (*XYSh)(1,J2,ISH);
  tauxip2 = xj2-xj;
  tauyip2 = yj2-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::onePointBackward(unsigned J, unsigned ISH)
{
  xj = (*XYSh)(0,J,ISH);
  yj = (*XYSh)(1,J,ISH);
  tauxim1 = xi-xj;
  tauyim1 = yi-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::twoPointsBackward(unsigned J2, unsigned ISH)
{
  xj2 = (*XYSh)(0,J2,ISH);
  yj2 = (*XYSh)(1,J2,ISH);
  tauxim2 = xj-xj2;
  tauyim2 = yj-yj2;
}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setTauIp1ToZero() {tauxip1 = 0; tauyip1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setTauIp2ToZero() {tauxip2 = 0; tauyip2 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setTauIm1ToZero() {tauxim1 = 0; tauyim1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4Pg::setTauIm2ToZero() {tauxim2 = 0; tauyim2 = 0;}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
