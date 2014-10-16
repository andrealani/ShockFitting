// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "RemeshingSF/CoNormPG.hh"
#include "RemeshingSF/ShpDpndnc.hh"
#include "Framework/Log.hh"

//----------------------------------------------------------------------------//

using namespace std;

//---------------------------------------------------------------------------//

namespace ShockFitting {

//---------------------------------------------------------------------------//

CoNormPG::CoNormPG(const std::string& objectName)
 :CoNorm("CoNormPG")
{
}

//----------------------------------------------------------------------------//

CoNormPG::~CoNormPG()
{
}

//----------------------------------------------------------------------------//

void CoNormPG::setup()
{
  LogToScreen(VERBOSE,"CoNormPG::setup() => start\n");

  LogToScreen(VERBOSE,"CoNormPG::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CoNormPG::unsetup()
{
  LogToScreen(VERBOSE,"CoNormPG::unsetup()\n");
}

//----------------------------------------------------------------------------//

void CoNormPG::remesh()
{
  LogToScreen(INFO,"CoNormPG::remesh()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setSize();

  // compute normal vector for each shock
  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {

    logfile("I: ",(*r_ZRoeShd)(0,I,ISH), " " ,(*r_ZRoeShd)(1,I,ISH));
    logfile("I: ",(*r_ZRoeShd)(2,I,ISH), " " ,(*r_ZRoeShd)(3,I,ISH));

     computeTau(ISH,I);

     // assign normal vector
     (*r_vShNor)(0,I,ISH) = -tauy;
     (*r_vShNor)(1,I,ISH) = taux;
   }
  }

  // compute normal vectors for typeSh="S"
  setVShNorForStype();

  // fix normal vectors for special points
  // it forces the direction of the contact discontinuity normal vector
  // in order to define an angle equal to 90° with mach stem
  // normal vector
  for (unsigned ISPPNTS=0; ISPPNTS<(*r_nSpecPoints); ISPPNTS++) {

   // special point: wall point without reflection
   if (r_typeSpecPoints->at(ISPPNTS)=="WPNRX") {setVShNorForWPNRX(ISPPNTS);}

   // special point: connection between two shocks
   else if (r_typeSpecPoints->at(ISPPNTS)=="C") {setVShNorForC(ISPPNTS);}

   // special point: triple point
   else if (r_typeSpecPoints->at(ISPPNTS)=="TP") {setVShNorForTP(ISPPNTS);}

   // write the computed normal vectors on tecplot file
   writeTecPlotFile();
  }

  logfile.Close();
}

//----------------------------------------------------------------------------//

void CoNormPG::computeTau(unsigned ISH, unsigned I)
{
  unsigned J, J2;
  ush = 0; vsh = 0;
  xi = (*r_XYSh)(0,I,ISH);
  yi = (*r_XYSh)(1,I,ISH);

  if (I < (r_nShockPoints->at(ISH)-1)) {
   // one point forward
   J=I+1;
   // coordinates of the one point forward
   onePointForward(J, ISH);

   // recover status for the forward point
   recoverStatus(J,ISH);

   if (I < (r_nShockPoints->at(ISH)-1)) {
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
  if (I!=0 && I!= r_nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepip(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depip1 = shockDepip.getDependence();
  }
  else if (I==r_nShockPoints->at(ISH)) {depip1=0; depim1=0;}

  if (I>2) {
   // one point backward
   J=I-1;
   // coordinates of the one point backward
   onePointBackward(J,ISH);

   // recover status for the backward point
   recoverStatus(J,ISH);

   if (I>3) {
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
  if (I!=0 && I!= r_nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepim(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depim1 = shockDepim.getDependence();
  }
  else if (I==0) {depip1=0; depim1=0;}

  setLp();
  setLm();

  if(I==0) {depim1=0; depip1=1; lm12=1;}
  if(I==r_nShockPoints->at(ISH)-1) {depim1=1; depip1=0; lp12=1.0;}

  taux = (depim1*tauxim1*lp12+depip1*tauxip1*lm12);
  tauy = (depim1*tauyim1*lp12+depip1*tauyip1*lm12);

  logfile(I, "tau ", taux, " ", tauy);
  logfile(" ", depim1, " ", depip1, "\n");

  tau = sqrt(taux*taux+tauy*tauy);
  taux = taux/tau;
  tauy = tauy/tau;
}

//----------------------------------------------------------------------------//

void CoNormPG::setVShNorForStype()
{
  for(unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   if(r_typeSh->at(ISH)=="S") {
    unsigned ii = 0;
    for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
     ui = (*r_ZRoeShd)(2,I,ISH)/(*r_ZRoeShd)(0,I,ISH);
     vi = (*r_ZRoeShd)(3,I,ISH)/(*r_ZRoeShd)(0,I,ISH);
     dum = ui*(*r_vShNor)(0,I,ISH)+vi*(*r_vShNor)(1,I,ISH);
     if(dum>0) {ii++;}
    }
    if (ii > r_nShockPoints->at(ISH)/2) {
     for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
      (*r_vShNor)(0,I,ISH) = -(*r_vShNor)(0,I,ISH);
      (*r_vShNor)(1,I,ISH) = -(*r_vShNor)(1,I,ISH);
     } // I
    } // if ii > r_nShockPoints->at(ISH)/2
   } // if r_typeSh->at(ISH)=="S"
  } // for
}

//----------------------------------------------------------------------------//

void CoNormPG::setVShNorForWPNRX(unsigned ISPPNTS)
{
  ISH1 = (*r_SHinSPPs)(0,0,ISPPNTS);
  I1 = (*r_SHinSPPs)(1,0,ISPPNTS)-1;
  IP1 = 1+I1*(r_nShockPoints->at(ISH1)-1);

  (*r_vShNor)(0,IP1,ISH1) = (*r_vShNor)(0,IP1,ISH1)/abs((*r_vShNor)(0,IP1,ISH1));
  (*r_vShNor)(1,IP1,ISH1) = 0;
}

//----------------------------------------------------------------------------//

void CoNormPG::setVShNorForC(unsigned ISPPNTS)
{
  ISH1 = (*r_SHinSPPs)(0,0,ISPPNTS);
  I1 = (*r_SHinSPPs)(1,0,ISPPNTS)-1;
  IP1 = 1+I1*(r_nShockPoints->at(ISH1)-1);

  ISH2 = (*r_SHinSPPs)(0,1,ISPPNTS);
  I2 = (*r_SHinSPPs)(1,1,ISPPNTS)-1;
  IP2 = 1+I2*(r_nShockPoints->at(ISH2)-1);

  nx1 = (*r_vShNor)(0,IP1,ISH1);
  nx1 = (*r_vShNor)(1,IP1,ISH1);
  nx2 = (*r_vShNor)(0,IP2,ISH2);
  nx2 = (*r_vShNor)(1,IP2,ISH2);

  nx1 = nx1+nx2;
  ny1 = ny1+ny2;

  dum = sqrt(nx1*nx1+ny1*ny1);
  nx1 = nx1/dum;
  ny1 = ny1/dum;

  (*r_vShNor)(0,IP1,ISH1) = nx1;
  (*r_vShNor)(1,IP1,ISH1) = ny1;
  (*r_vShNor)(0,IP2,ISH1) = nx1;
  (*r_vShNor)(1,IP2,ISH1) = ny1;
}

//----------------------------------------------------------------------------//

void CoNormPG::setVShNorForTP(unsigned ISPPNTS)
{
  // define shocks and edge indeces

  // incident shock
  ISH1 = (*r_SHinSPPs)(0,0,ISPPNTS);
  I1 = (*r_SHinSPPs)(1,0,ISPPNTS)-1;
  IP1 = 1+I1*(r_nShockPoints->at(ISH1)-1);
  // reflected shock
  ISH2 = (*r_SHinSPPs)(0,1,ISPPNTS);
  I2 = (*r_SHinSPPs)(1,1,ISPPNTS)-1;
  IP2 = 1+I2*(r_nShockPoints->at(ISH2)-1);
  // Mach stem
  ISH3 = (*r_SHinSPPs)(0,2,ISPPNTS);
  I3 = (*r_SHinSPPs)(1,2,ISPPNTS)-1;
  IP3 = 1+I3*(r_nShockPoints->at(ISH3)-1);
  // contact discontinuity
  ISH4 = (*r_SHinSPPs)(0,3,ISPPNTS);
  I4 = (*r_SHinSPPs)(1,3,ISPPNTS)-1;
  IP4 = 1+I4*(r_nShockPoints->at(ISH4)-1);

  nx2 = (*r_vShNor)(0,IP2,ISH2);
  ny2 = (*r_vShNor)(1,IP2,ISH2);

  nx4 = (*r_vShNor)(0,IP4,ISH4);
  ny4 = (*r_vShNor)(1,IP4,ISH4);

  dum = nx2*nx4+ny2*ny4;
  if (dum < 0) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH4); I++) {
    (*r_vShNor)(0,I,ISH4) = -(*r_vShNor)(0,I,ISH4);
    (*r_vShNor)(1,I,ISH4) = -(*r_vShNor)(1,I,ISH4);
   }
  }
}

//----------------------------------------------------------------------------//

void CoNormPG::writeTecPlotFile()
{
  ofstream tecfile;
  tecfile.open("shocknor.dat");

  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   tecfile << "TITLE = Shock normals\n";
   tecfile << "VARIABLES = X Y Z(1) Z(2) NX NY\n";
   tecfile << "ZONE T='sampletext' F=FEPOINT ET=TRIANGLE ";
   tecfile << "N= " << r_nShockPoints->at(ISH);
   tecfile << "E= " << r_nShockPoints->at(ISH)-1 << "\n";
   for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
    for (unsigned K=0; K<(*ndim); K++) {tecfile << (*r_XYSh)(K,I,ISH) << " ";}
    tecfile << 1 << " " << 1 << " ";
    for (unsigned K=0; K<(*ndim); K++) {tecfile << (*r_vShNor)(K,I,ISH) << " ";}
   }
   for (unsigned I=0; I<r_nShockPoints->at(ISH)-1; I++) {
    tecfile << I << " " << I+1 << " " << I << "\n";
   }
  }
}

//----------------------------------------------------------------------------//

void CoNormPG::recoverStatus(unsigned J, unsigned ISH)
{
  uj = (*r_ZRoeShd)(2,J,ISH)/(*r_ZRoeShd)(0,J,ISH);
  vj = (*r_ZRoeShd)(3,J,ISH)/(*r_ZRoeShd)(0,J,ISH);
  roj = (*r_ZRoeShd)(0,J,ISH)*(*r_ZRoeShd)(0,J,ISH);
  help = ((*r_ZRoeShd)(2,J,ISH),2)+pow((*r_ZRoeShd)(3,J,ISH),2);
  pj = (*gm1)/(*gam) * ((*r_ZRoeShd)(0,J,ISH)*(*r_ZRoeShd)(1,J,ISH)-0.5*help);
  aj = sqrt((*gam)*pj/roj);

  if (r_typeSh->at(ISH)=="D") { aj = 0;}
}

//----------------------------------------------------------------------------//

void CoNormPG::setLp()
{
  lp12 = pow(tauxip1,2)+pow(tauyip1,2);
  lp22 = pow(tauxip2,2)+pow(tauyip2,2);
  lp1 = sqrt(lp12);
  lp2 = sqrt(lp22); 
}

//----------------------------------------------------------------------------//

void CoNormPG::setLm()
{
  lm12 = pow(tauxim1,2)+pow(tauyim1,2);
  lm22 = pow(tauxim2,2)+pow(tauyim2,2);
  lm1 = sqrt(lm12);
  lm2 = sqrt(lp22);
}

//----------------------------------------------------------------------------//

void CoNormPG::onePointForward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxip1 = xj-xi;
  tauyip1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNormPG::twoPointsForward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxip2 = xj2-xi;
  tauyip2 = yj2-yi;
}

//----------------------------------------------------------------------------//

void CoNormPG::onePointBackward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxim1 = xj-xi;
  tauyim1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNormPG::twoPointsBackward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxim2 = xj2-xi;
  tauyim2 = yj2-yi;
}

//----------------------------------------------------------------------------//

void CoNormPG::setTauIp1ToZero() {tauxip1 = 0; tauyip1 = 0;}

//----------------------------------------------------------------------------//

void CoNormPG::setTauIp2ToZero() {tauxip2 = 0; tauyip2 = 0;}

//----------------------------------------------------------------------------//

void CoNormPG::setTauIm1ToZero() {tauxim1 = 0; tauyim1 = 0;}

//----------------------------------------------------------------------------//

void CoNormPG::setTauIm2ToZero() {tauxim2 = 0; tauyim2 = 0;}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
