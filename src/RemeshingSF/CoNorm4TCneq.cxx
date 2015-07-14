// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/CoNorm4TCneq.hh"
#include "RemeshingSF/ShpDpndnc.hh"
#include "Framework/MeshData.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//---------------------------------------------------------------------------//

namespace ShockFitting {

//---------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CoNorm4TCneq, Remeshing> 
computeNormalVector4TCneqProv("CoNorm4TCneq");

//--------------------------------------------------------------------------//

CoNorm4TCneq::CoNorm4TCneq(const std::string& objectName) :
  Remeshing(objectName)
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
  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    unsigned shpoin = I+1; // c++ indeces start from 0
    logfile(shpoin, " ");
    for(unsigned J=0; J<(*ndof); J++) {
     logfile((*ZRoeShd)(J,I,ISH), " ");
    } // J
    logfile("\n");
   } // I
   logfile("\n");
  } // ISH
  logfile("\n");

  // compute normal vector for each shock
  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {

    computeTau(ISH,I);

    // assign normal vector
    (*vShNor)(0,I,ISH)=tauy;
    (*vShNor)(1,I,ISH)=-taux;
   }
  }


  // compute normal vectors for typeSh="S"
  setVShNorForStype();

  // fix normal vectors for special points
  // it forces the direction of the contact discontinuity normal vector
  // in order to define an angle equal to 90° with mach stem normal vectr
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

  // de-allocate ZRoeShd
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::computeTau(unsigned ISH, unsigned I)
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
   ipoin = J+1; // ipoin is used in the log file (c++ indeces start from 0)
   recoverState("forward",I,J,ISH);

   if (I < (nShockPoints->at(ISH)-2)) {
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
//  else if (I==nShockPoints->at(ISH)-1) {depip1=0; depim1=1;}


  if (I>0) {
   // one point backward
   J=I-1;
   // coordinates of the one point backward
   onePointBackward(J,ISH);

   // recover status for the backward point
   recoverState("backward",I,J,ISH);

   ipoin = J+1; // ipoin is used in the log file (c++ indeces start from 0)

   if (I>1) {
    // two points backward
    J2=I-2;
    // coordinates of two points backward
    twoPointsBackward(J2,ISH);
   }
   else {setTauIm2ToZero();}
  }
  else { // if I=0
   setTauIm1ToZero();
   setTauIm2ToZero();
  }
  
  // evaluate dependency
  // if I is the first point it is computed in forward direction
  if (I!=0 && I!= nShockPoints->at(ISH)-1) {
   ShpDpndnc shockDepim(xi,yi,ush,vsh,xj,yj,uj,vj,aj);
   depim1 = shockDepim.getDependence();
  }
//  else if (I==0) {depip1=1; depim1=0;}

  setLp();
  setLm();

  if(I==0) {depim1=0; depip1=1; lm12=1.0;}
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

void CoNorm4TCneq::setVShNorForStype()
{
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   if(typeSh->at(ISH)=="S") {
    unsigned ii = 0;
    for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
     unsigned zrho = 0;
     for(unsigned ISP=0; ISP<(*nsp); ISP++) {
      zrho = zrho + (*ZRoeShd)(ISP,I,ISH);
     }
     ui = (*ZRoeShd)((*IX),I,ISH)/zrho;
     vi = (*ZRoeShd)((*IY),I,ISH)/zrho;
     dum = ui * (*vShNor)(0,I,ISH) + vi * (*vShNor)(1,I,ISH);
     if (dum>0) {++ii;}
    }
    if (ii > nShockPoints->at(ISH)/2) {
     for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
      (*vShNor)(0,I,ISH) = -(*vShNor)(0,I,ISH);
      (*vShNor)(1,I,ISH) = -(*vShNor)(1,I,ISH);
     }
    }
   } // I
  } // for
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  (*vShNor)(0,IP.at(0),ISH.at(0)) =
      (*vShNor)(0,IP.at(0),ISH.at(0))/abs((*vShNor)(0,IP.at(0),ISH.at(0)));
  (*vShNor)(1,IP.at(0),ISH.at(0)) = 0;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForC(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
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
  (*vShNor)(0,IP.at(1),ISH.at(1)) = nx1;
  (*vShNor)(1,IP.at(1),ISH.at(1)) = ny1;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setVShNorForTP(unsigned ISPPNTS)
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

void CoNorm4TCneq::writeTecPlotFile()
{
  FILE* tecfile;
  tecfile = fopen("shocknor.dat","w");

  for (unsigned ISH=0; ISH<(*nShocks); ISH++) {
   fprintf(tecfile,"%s","TITLE = Shock normals\n");
   fprintf(tecfile,"%s","VARIABLES = X Y Z(1) Z(2) NX NY\n");
   fprintf(tecfile,"%s","ZONE T='sampletext', F = FEPOINT, ET = TRIANGLE ");
   fprintf(tecfile,"%s %5u","N = ",nShockPoints->at(ISH));
   fprintf(tecfile,"%s %5u %s",", E = ",nShockPoints->at(ISH)-1,"\n");
   for (unsigned I=0; I<nShockPoints->at(ISH); I++) {
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++)
     {fprintf(tecfile,"%32.16E %s",(*XYSh)(K,I,ISH)," ");}
    fprintf(tecfile,"%s","\n");
    fprintf(tecfile,"%s","1  1 ");
    for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++)
     {fprintf(tecfile,"%32.16E %s",(*vShNor)(K,I,ISH)," ");}
    fprintf(tecfile,"%s","\n");
   }
   for (unsigned I=0; I<nShockPoints->at(ISH)-1; I++) {
    fprintf(tecfile,"%u %s %u %s %u %s",I+1," ",I+2," ",I+1,"\n");
   }
  }

  fclose(tecfile);
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::recoverState(string direction, unsigned I,
				 unsigned J, unsigned ISH)
{
  double zrho = 0;
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   zrho = zrho + (*ZRoeShd)(ISP,J,ISH);
  }
  roj = zrho*zrho;
  double rhoHf = 0; double Rg = 0; double Cv = 0;
  for (unsigned ISP=0; ISP<(*nsp); ISP++) {
   if(direction=="forward") { rhoHf = rhoHf + (*ZRoeShd)(ISP,J,ISH) * hf->at(ISP); }
   if(direction=="backward") { rhoHf = (*ZRoeShd)(ISP,J,ISH) * hf->at(ISP); }
   Rg = Rg + (*ZRoeShd)(ISP,J,ISH)*Rs->at(ISP);
   Cv = Cv + (*ZRoeShd)(ISP,J,ISH)*Rs->at(ISP)/(gams->at(ISP)-1);
  }
  rhoHf = rhoHf*zrho;
  Rg = Rg*zrho;
  Cv = Cv*zrho;

  double gamm1, gammam, gm1oga;
  gamm1 = Rg/Cv;
  gammam = 1+ gamm1;
  gm1oga = gamm1/gammam;

  uj = (*ZRoeShd)((*IX),J,ISH)/zrho;
  vj = (*ZRoeShd)((*IY),J,ISH)/zrho;
  help = pow((*ZRoeShd)((*IX),J,ISH),2)+pow((*ZRoeShd)((*IY),J,ISH),2);
  pj = gm1oga * ( zrho * ((*ZRoeShd)((*IE),J,ISH)) -
                  0.5e0 * help - rhoHf -
                  zrho * ((*ZRoeShd)((*IEV),J,ISH)));
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

  if (typeSh->at(ISH)=="D") { aj = 0;}
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
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
  xj = (*XYSh)(0,J,ISH);
  yj = (*XYSh)(1,J,ISH);
  tauxip1 = xj-xi;
  tauyip1 = yj-yi;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::twoPointsForward(unsigned J2, unsigned ISH)
{
  xj2 = (*XYSh)(0,J2,ISH);
  yj2 = (*XYSh)(1,J2,ISH);
  tauxip2 = xj2-xj;
  tauyip2 = yj2-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::onePointBackward(unsigned J, unsigned ISH)
{
  xj = (*XYSh)(0,J,ISH);
  yj = (*XYSh)(1,J,ISH);
  tauxim1 = xi-xj;
  tauyim1 = yi-yj;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::twoPointsBackward(unsigned J2, unsigned ISH)
{
  xj2 = (*XYSh)(0,J2,ISH);
  yj2 = (*XYSh)(1,J2,ISH);
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

void CoNorm4TCneq::setAddress()
{
  unsigned start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
                   PhysicsInfo::getnbShPointsMax() *
                   PhysicsInfo::getnbShMax() *
                   PhysicsInfo::getnbDofMax();
  ZRoeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroe->at(start));
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setSize()
{
  vShNor->resize(PhysicsInfo::getnbDim(),
                 PhysicsInfo::getnbShPointsMax(),
                 PhysicsInfo::getnbShMax());
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::freeArray()
{
  delete ZRoeShd;
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroe = MeshData::getInstance().getData <vector<double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void CoNorm4TCneq::setPhysicsData()
{
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  Rs = PhysicsData::getInstance().getData <vector<double> > ("RS");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
      PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nSpecPoints =
      PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  typeSh =
      PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  typeSpecPoints =
      PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  XYSh = PhysicsData::getInstance().getData <Array3D<double> > ("XYSH");
  vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  SHinSPPs =
      PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
