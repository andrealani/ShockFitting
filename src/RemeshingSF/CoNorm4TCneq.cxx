// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "RemeshingSF/CoNorm4TCneq.hh"
#include "RemeshingSF/ShpDpndnc.hh"
#include "Framework/Log.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//---------------------------------------------------------------------------//

namespace ShockFitting {

//---------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CoNorm4TCneq, CoNorm> 
computeNormalVector4TCneqProv("CoNorm4TCneq");

//--------------------------------------------------------------------------//

CoNorm4TCneq::CoNorm4TCneq(const std::string& objectName) :
  CoNorm(objectName)
{
}

//----------------------------------------------------------------------------//

CoNorm4TCneq::~CoNorm4TCneq()
{
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setup()
{
  LogToScreen(VERBOSE,"CoNorm4TCneq::setup() => start\n");

  LogToScreen(VERBOSE,"CoNorm4TCneq::setup() => end\n");
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::unsetup()
{
  LogToScreen(VERBOSE,"CoNorm4TCneq::unsetup()\n");
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::remesh()
{
  LogToScreen(INFO,"CoNorm4TCneq::remesh()\n");

  logfile.Open(getClassName());

  setMeshData();
  setPhysicsData();

  setAddress();

  setSize();

  // write downstream status on log file
  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
    unsigned shpoin = I+1; // c++ indeces start from 0
    logfile(shpoin, " ");
    for(unsigned J=0; J<(*ndof); J++) {
     logfile((*r_ZRoeShd)(J,I,ISH), " ");
    } // J
    logfile("\n");
   } // I
   logfile("\n");
  } // ISH
  logfile("\n");

  // compute normal vector for each shock
  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {

    computeTau(ISH,I);

    // assign normal vector
    (*r_vShNor)(0,I,ISH)=tauy;
    (*r_vShNor)(1,I,ISH)=-taux;
   }
  }

  // compute normal vectors for typeSh="S"
  setVShNorForStype();

  // fix normal vectors for special points
  // it forces the direction of the contact discontinuity normal vector
  // in order to define an angle equal to 90° with mach stem normal vectr
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

void CoNorm4TCneq::computeTau(unsigned ISH, unsigned I)
{
  unsigned J, J2, ipoin;
  ush = 0; vsh = 0;
  xi = (*r_XYSh)(0,I,ISH);
  yi = (*r_XYSh)(1,I,ISH);

  if (I < (r_nShockPoints->at(ISH)-1)) {
   // one point forward
   J=I+1;
   // coordinates of the one point forward
   onePointForward(J, ISH);

   // recover status for the forward point
   ipoin = J+1; // c++ indeces start from 0
   recoverStatus("forward",I,J,ISH);

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
  else if (I==r_nShockPoints->at(ISH)) {depip1=0; depim1=1;}


  if (I>0) {
   // one point backward
   J=I-1;
   // coordinates of the one point backward
   onePointBackward(J,ISH);

   // recover status for the backward point
   ipoin = J+1; // c++ indeces start from 0
   recoverStatus("backward",I,J,ISH);

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
  if (I!=0 && I!= r_nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepim(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depim1 = shockDepim.getDependence();
  }
  else if (I==0) {depip1=1; depim1=0;}

  setLp();
  setLm();

  if(I==0) {depim1=0; depip1=1; lm12=1;}
  if(I==r_nShockPoints->at(ISH)-1) {depim1=1; depip1=0; lp12=1.0;}

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

void CoNorm4TCneq::setVShNorForStype()
{
  for(unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   if(r_typeSh->at(ISH)=="S") {
    unsigned ii = 0;
    for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
     unsigned zrho = 0;
     for(unsigned ISP=0; ISP<(*nsp); ISP++) {
      zrho = zrho + (*r_ZRoeShd)(ISP,I,ISH);
     }
     ui = (*r_ZRoeShd)((*ix),I,ISH)/zrho;
     vi = (*r_ZRoeShd)((*iy),I,ISH)/zrho;
     dum = ui * (*r_vShNor)(0,I,ISH) + vi * (*r_vShNor)(1,I,ISH);
     if (dum>0) {++ii;}
    }
    if (ii < r_nShockPoints->at(ISH)/2 - 1) { break; }
     for (unsigned I=0; I<r_nShockPoints->at(ISH); I++) {
      (*r_vShNor)(0,I,ISH) = -(*r_vShNor)(0,I,ISH);
      (*r_vShNor)(1,I,ISH) = -(*r_vShNor)(1,I,ISH);
     } // I
   } // if r_typeSh->at(ISH)=="S"
  } // for
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  (*r_vShNor)(0,IP.at(0),ISH.at(0)) =
      (*r_vShNor)(0,IP.at(0),ISH.at(0))/abs((*r_vShNor)(0,IP.at(0),ISH.at(0)));
  (*r_vShNor)(1,IP.at(0),ISH.at(0)) = 0;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForC(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  nx1 = (*r_vShNor)(0,IP.at(0),ISH.at(0));
  ny1 = (*r_vShNor)(1,IP.at(0),ISH.at(0));
  nx2 = (*r_vShNor)(0,IP.at(1),ISH.at(1));
  ny2 = (*r_vShNor)(1,IP.at(1),ISH.at(1));

  nx1 = nx1+nx2;
  ny1 = ny1+ny2;

  dum = sqrt(nx1*nx1+ny1*ny1);
  nx1 = nx1/dum;
  ny1 = ny1/dum;

  (*r_vShNor)(0,IP.at(0),ISH.at(0)) = nx1;
  (*r_vShNor)(1,IP.at(0),ISH.at(0)) = ny1;
  (*r_vShNor)(0,IP.at(1),ISH.at(1)) = nx1;
  (*r_vShNor)(1,IP.at(1),ISH.at(1)) = ny1;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForTP(unsigned ISPPNTS)
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

void CoNorm4TCneq::writeTecPlotFile()
{
  ofstream tecfile;
  tecfile.open("shocknor.dat");

  for (unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   tecfile << "TITLE = Shock normals\n";
   tecfile << "VARIABLES = X Y Z(1) Z(2) NX NY\n";
   tecfile << "ZONE T='sampletext', F = FEPOINT, ET = TRIANGLE ";
   tecfile << "N = " << r_nShockPoints->at(ISH);
   tecfile << ", E = " << r_nShockPoints->at(ISH)-1 << "\n";
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

void CoNorm4TCneq::recoverStatus(string direction, unsigned I,
				 unsigned J, unsigned ISH)
{
  double zrho = 0;
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   zrho = zrho + (*r_ZRoeShd)(ISP,J,ISH);
  }
  roj = zrho*zrho;

  double rhoHf = 0; double Rg = 0; double Cv = 0;
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   rhoHf = rhoHf + (*r_ZRoeShd)(ISP,J,ISH) * hf->at(ISP);
   Rg = Rg + (*r_ZRoeShd)(ISP,J,ISH)*Rs->at(ISP);
   Cv = Cv + (*r_ZRoeShd)(ISP,J,ISH)*Rs->at(ISP)/(gams->at(ISP)-1);
  }
  rhoHf = rhoHf*zrho;
  Rg = Rg*zrho;
  Cv = Cv*zrho;

  double gamm1, gammam, gm1oga;
  gamm1 = Rg/Cv;
  gammam = 1+ gamm1;
  gm1oga = gamm1/gammam;

  uj = (*r_ZRoeShd)((*ix),J,ISH)/zrho;
  vj = (*r_ZRoeShd)((*iy),J,ISH)/zrho;
  help = pow((*r_ZRoeShd)((*ix),J,ISH),2)+pow((*r_ZRoeShd)((*iy),J,ISH),2);
  pj = gm1oga* (zrho*(*r_ZRoeShd)((*ie),J,ISH)-
       0.5 * help - rhoHf - zrho * (*r_ZRoeShd)((*iev),J,ISH));
  aj = sqrt(gammam*pj/roj);

  if(direction=="forward") {  
   logfile("g = ",gammam,"\n");
   logfile("(g-1)/g = ", gm1oga, "\n");
   logfile("rho = ", roj, "\n");
   logfile("gref = ", (*gref),"\n");
   logfile("RHOHF = ",rhoHf, "\n");
   logfile("RG = ", Rg, "\n");
   logfile("Cv = ",Cv, "\n");

   unsigned ipoin = I+1; //c++ indeces start from 0
   unsigned ish = ISH+1; //c++ indeces start from 0
   logfile ("\nShock no. ", ish, "point no. ", ipoin, " ");
   logfile (pj, " ", aj , "\n");
  }

  if (r_typeSh->at(ISH)=="D") { aj = 0;}
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
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

void CoNorm4TCneq::setLp()
{
  lp12 = pow(tauxip1,2)+pow(tauyip1,2);
  lp22 = pow(tauxip2,2)+pow(tauyip2,2);
  lp1 = sqrt(lp12);
  lp2 = sqrt(lp22);
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setLm()
{
  lm12 = pow(tauxim1,2)+pow(tauyim1,2);
  lm22 = pow(tauxim2,2)+pow(tauyim2,2);
  lm1 = sqrt(lm12);
  lm2 = sqrt(lp22);
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::onePointForward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxip1 = xj-xi;
  tauyip1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::twoPointsForward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxip2 = xj2-xj;
  tauyip2 = yj2-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::onePointBackward(unsigned J, unsigned ISH)
{
  xj = (*r_XYSh)(0,J,ISH);
  yj = (*r_XYSh)(1,J,ISH);
  tauxim1 = xi-xj;
  tauyim1 = yi-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::twoPointsBackward(unsigned J2, unsigned ISH)
{
  xj2 = (*r_XYSh)(0,J2,ISH);
  yj2 = (*r_XYSh)(1,J2,ISH);
  tauxim2 = xj-xj2;
  tauyim2 = yj-yj2;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setTauIp1ToZero() {tauxip1 = 0; tauyip1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setTauIp2ToZero() {tauxip2 = 0; tauyip2 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setTauIm1ToZero() {tauxim1 = 0; tauyim1 = 0;}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setTauIm2ToZero() {tauxim2 = 0; tauyim2 = 0;}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
