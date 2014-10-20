// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "RemeshingSF/CoNorm4B.hh"
#include "Framework/Log.hh"

//----------------------------------------------------------------------------//

using namespace std;

//---------------------------------------------------------------------------//

namespace ShockFitting {

//---------------------------------------------------------------------------//

CoNorm4B::CoNorm4B(const std::string& objectName)
 :CoNorm("CoNorm4B")
{
}

//----------------------------------------------------------------------------//

CoNorm4B::~CoNorm4B()
{
}

//----------------------------------------------------------------------------//

void CoNorm4B::setup()
{
  LogToScreen(VERBOSE,"CoNorm4B::setup() => start\n");

  LogToScreen(VERBOSE,"CoNorm4B::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CoNorm4B::unsetup()
{
  LogToScreen(VERBOSE,"CoNorm4B::unsetup()\n");
}

//----------------------------------------------------------------------------//

void CoNorm4B::remesh()
{
  LogToScreen(INFO,"CoNorm4B::remesh()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setSize();

  // compute normal vector for each shock
  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {

    logfile("Shock variables ",(*r_ZRoeShd)(0,I,ISH));

    computeTau(ISH,I);

    // assign normal vector
    (*r_vShNor)(0,I,ISH) = tauy;
    (*r_vShNor)(1,I,ISH) = -taux;
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

void CoNorm4B::computeTau(unsigned ISH,unsigned I)
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
   if (I < (r_nShockPoints->at(ISH)-1)) {
    // two points forward
    J2=I+2;
    // coordinates of the two points forward
    twoPointsForward(J2,ISH);
   } // J2
   else { setTauIp2ToZero();}
  } // J
  else { // I = nShockPoints(ISH)-1 
   setTauIp1ToZero();
   setTauIp2ToZero();
  }

  if (I>0) {
   // one point backward
   J=I-1;
   // coordinates of the one point backward
   onePointBackward(J,ISH);
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

  setLp();
  setLm();

  taux = tauxim1*(lm12+lm22)-tauxim2*lm12;
  tauy = tauyim1*(lm12+lm22)-tauyim2*lm12;
  if( I == 0 ){taux = tauxip1; tauy = tauyip1;}

  tau = sqrt(taux*taux+tauy*tauy);
  taux = taux/tau;
  tauy = tauy/tau;
 
  logfile(I, " ", I, ": ");
  logfile(tau, " ", taux, " ", tauy, "\n");  
  logfile(I, " ", J, ": ");
  logfile(lm12, " ", tauxim1, " ", tauyim1, "\n");
  logfile(I, " ", J2, ": ");
  logfile(lm22, " ", tauxim2, " ", tauyim2, "\n");
}

//----------------------------------------------------------------------------//

void CoNorm4B::setVShNorForStype()
{
  for(unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   if(r_typeSh->at(ISH)=="S") {
    unsigned ii = 0;
    for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {++ii;}
    if (ii > r_nShockPoints->at(ISH)/2) {
     for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
      (*r_vShNor)(0,I,ISH) = -(*r_vShNor)(0,I,ISH);
      (*r_vShNor)(1,I,ISH) = -(*r_vShNor)(1,I,ISH);
     } // I
    } // if ii > r_nShockPoints->at(ISH)/2
   } // if r_typeSh->at(ISH)=="S"
  } // for ISH
}

//----------------------------------------------------------------------------//

void CoNorm4B::setVShNorForWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  (*r_vShNor)(0,IP.at(0),ISH.at(0)) =
      (*r_vShNor)(0,IP.at(0),ISH.at(0))/abs((*r_vShNor)(0,IP.at(0),ISH.at(0)));
  (*r_vShNor)(1,IP.at(0),ISH.at(0)) = 0;
}

//----------------------------------------------------------------------------//

void CoNorm4B::setVShNorForC(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  nx1 = (*r_vShNor)(0,IP.at(0),ISH.at(0));
  nx1 = (*r_vShNor)(1,IP.at(0),ISH.at(0));
  nx2 = (*r_vShNor)(0,IP.at(1),ISH.at(1));
  ny2 = (*r_vShNor)(1,IP.at(1),ISH.at(1));

  nx1 = nx1+nx2;
  ny1 = ny1+ny2;

  dum = sqrt(nx1*nx1+ny1*ny1);
  nx1 = nx1/dum;
  ny1 = ny1/dum;

  (*r_vShNor)(0,IP.at(0),ISH.at(0)) = nx1;
  (*r_vShNor)(1,IP.at(0),ISH.at(0)) = ny1;
  (*r_vShNor)(0,IP.at(1),ISH.at(0)) = nx1;
  (*r_vShNor)(1,IP.at(1),ISH.at(0)) = ny1;

}

//----------------------------------------------------------------------------//

void CoNorm4B::setVShNorForTP(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock 
  // ISH.at(1) reflected shock
  // ISH.at(2) Mach stem
  // ISH.at(3) contact discontinuity
  setShockIndeces(4,ISPPNTS);

  nx2 = (*r_vShNor)(0,IP.at(1),ISH.at(1));
  ny2 = (*r_vShNor)(1,IP.at(1),ISH.at(1));

  nx4 = (*r_vShNor)(0,IP.at(3),ISH.at(3));
  ny4 = (*r_vShNor)(1,IP.at(3),ISH.at(3));

  dum = nx2*nx4+ny2*ny4;
  if (dum < 0) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH.at(3)); I++) {
    (*r_vShNor)(0,I,ISH.at(3)) = -(*r_vShNor)(0,I,ISH.at(3));
    (*r_vShNor)(1,I,ISH.at(3)) = -(*r_vShNor)(1,I,ISH.at(3));
   }
  }
}

//----------------------------------------------------------------------------//

void CoNorm4B::writeTecPlotFile()
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

void CoNorm4B::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
{
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*r_SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*r_SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (r_nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//----------------------------------------------------------------------------//

void CoNorm4B::onePointForward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxip1 = xj-xi;
  tauyip1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNorm4B::twoPointsForward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxip2 = xj2-xj;
  tauyip2 = yj2-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4B::onePointBackward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxim1 = xi-xj;
  tauyim1 = yi-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4B::twoPointsBackward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxim2 = xj-xj2;
  tauyim2 = yj-yj2;
}

//----------------------------------------------------------------------------//

void CoNorm4B::setLm()
{
  lm12=pow(tauxim1,2)+pow(tauyim1,2);
  lm22=pow(tauxim2,2)+pow(tauyim2,2);

}

//----------------------------------------------------------------------------//

void CoNorm4B::setLp()
{
  lp12=pow(tauxip1,2)+pow(tauyip1,2);
  lp22=pow(tauxip2,2)+pow(tauyip2,2);
}

//----------------------------------------------------------------------------//

void CoNorm4B::setTauIp1ToZero() {tauxip1 = 0; tauyip1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4B::setTauIp2ToZero() {tauxip2 = 0; tauyip2 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4B::setTauIm1ToZero() {tauxim1 = 0; tauyim1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4B::setTauIm2ToZero() {tauxim2 = 0; tauyim2 = 0;}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
