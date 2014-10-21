// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/FxMshSps.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "RemeshingSF/FindBEdg.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<FxMshSps, Remeshing> fixMshSpsProv("FxMshSps");

//--------------------------------------------------------------------------//

FxMshSps::FxMshSps (const std::string& objectName)
 : Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FxMshSps::~FxMshSps()
{
}

//--------------------------------------------------------------------------//

void FxMshSps::setup()
{
  LogToScreen(VERBOSE,"FxMshSps::setup() => start\n");

  logfile.Open("FxMshSps");

  LogToScreen(VERBOSE,"FxMshSps::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FxMshSps::unsetup()
{
  logfile.Close();

  LogToScreen(VERBOSE,"FxMshSps::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FxMshSps::remesh()
{
  LogToScreen(INFO, "FxMshSps::remesh()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  createShockNodesList();

  /// create downstream shock edges
  createDownShEdges();

  for(unsigned ISPPNTS=0; ISPPNTS<(*r_nSpecPoints); ISPPNTS++) {

   if       (r_typeSpecPoints->at(ISPPNTS) == "IPX"   ||
             r_typeSpecPoints->at(ISPPNTS) == "IPY"   ||
             r_typeSpecPoints->at(ISPPNTS) == "OPX"   ||
             r_typeSpecPoints->at(ISPPNTS) == "OPY"   ||
             r_typeSpecPoints->at(ISPPNTS) == "WPNRX" ||
             r_typeSpecPoints->at(ISPPNTS) == "WPNRY") 
      { fixMeshForIPorOPorWPNR(ISPPNTS); }

    else if (r_typeSpecPoints->at(ISPPNTS) == "RRX")
      { fixMeshForRRX(ISPPNTS); }

    else if (r_typeSpecPoints->at(ISPPNTS) == "TP" ||
	     r_typeSpecPoints->at(ISPPNTS) == "QP" ||
             r_typeSpecPoints->at(ISPPNTS) == "EP" ||
             r_typeSpecPoints->at(ISPPNTS) == "C" )
       {} // nothing to do
    else 
       { cout << "Condition not implemented\n"; exit(1);}
  }

  (*nbfacSh) = ibfac;
}

//--------------------------------------------------------------------------//

void FxMshSps::createShockNodesList()
{
  ISHPlistu.resize((*npshmax),(*nshmax));
  ISHPlistd.resize((*npshmax),(*nshmax));
  unsigned ilist = *npoin;
  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) {
    ++ilist;
    ISHPlistu(I,ISH) = ilist;
   }
  }

  for(unsigned ISH=0; ISH<(*nshmax); ISH++) {
   for(unsigned I=0; I<(*npshmax); I++) { 
    ++ilist;
    ISHPlistd(I,ISH) = ilist;
   }
  }
}

//--------------------------------------------------------------------------//

void FxMshSps::createDownShEdges()
{
  ibfac = (*nbfac);
  for(unsigned ISH=0; ISH<(*r_nShocks); ISH++) {
   r_nShockEdges->at(ISH) = r_nShockPoints->at(ISH)-1;
   for (unsigned I=0; I<(r_nShockEdges->at(ISH)); I++) {
    (*bndfac)(0,ibfac) = ISHPlistd(I,ISH);
    (*bndfac)(1,ibfac) = ISHPlistd(I+1,ISH);
    (*bndfac)(2,ibfac) = 10;
    ++ibfac;
   }

   for (unsigned I=0; I<(r_nShockEdges->at(ISH)); I++) {
    (*bndfac)(0,ibfac) = ISHPlistu(I,ISH);
    (*bndfac)(1,ibfac) = ISHPlistu(I+1,ISH);
    (*bndfac)(2,ibfac) = 10;
    ++ibfac;
   }
  }
}

//--------------------------------------------------------------------------//

void FxMshSps::fixMeshForIPorOPorWPNR(unsigned ISPPNTS)
{

  setShockIndeces(1,ISPPNTS);
  unsigned ishock = ISH.at(0)+1 ;

  FindBEdg lookForBE;

  // look for boundary edge for downstream point
  setShockPointCoor("Down",IP.at(0),ISH.at(0));

  logfile("\nTypeSpecPoints: ", r_typeSpecPoints->at(ISPPNTS), "\n");

  iedg1 = lookForBE.getBEdg(xsh, ysh);
  s1 = lookForBE.getS();

  if (iedg1==-1) {
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  int iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(1): ", s1, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg1), " ", (*bndfac)(1,iedg1),"\n");

  // look for boundary edge for upstream point
  setShockPointCoor("Up",IP.at(0),ISH.at(0)); 
  iedg2 = lookForBE.getBEdg(xsh, ysh);
  s2 = lookForBE.getS();

  if (iedg2==-1) { 
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  iedg = iedg2+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(2): ", s2, " ");
  logfile(xsh," ", ysh, " ", iedg); 
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg2), " ", (*bndfac)(1,iedg2),"\n");

  // check if the shock crosses a boundary point
  checkShockBndryEdgeCrossing();

  iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile ("****************************");
  logfile ("    Shock: ", ishock, "\n");
  logfile ("    ", iedg, "\n");
  logfile ("****************************");


  // Split the existing edges of the background mesh
  // and create 2 new edges at  boundary
  // if the first point of the shock  is at the boundary
  if (iedg1>0) { splitEdges();
		 createNewEdges(IP.at(0), ISH.at(0)); }
}

//--------------------------------------------------------------------------//

void FxMshSps::fixMeshForRRX(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);
  unsigned ishock = ISH.at(1)+1;

  FindBEdg lookForBE;

  // look for boundary edge for downstream point
  setShockPointCoor("Down",IP.at(1),ISH.at(1));
  iedg1 = lookForBE.getBEdg(xsh, ysh);
  s1 = lookForBE.getS();
  if (iedg1==-1) { 
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  int iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(1): ", s1, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg1), " ", (*bndfac)(1,iedg1),"\n");

  // look for boundary edge for upstream point
  setShockPointCoor("Up",IP.at(0),ISH.at(0));
  iedg2 = lookForBE.getBEdg(xsh, ysh);
  s2 = lookForBE.getS();
  if (iedg2==-1) { 
   ishock = ISH.at(0)+1;
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  iedg = iedg2+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(2): ", s2, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg2), " ", (*bndfac)(1,iedg2),"\n");

  // check if the shock crosses a boundary point
  checkShockBndryEdgeCrossing();

  iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile ("****************************");
  logfile ("    Shock: ", ishock, "\n");
  logfile ("    ", iedg, "\n");
  logfile ("****************************");

  // Split the existing edges of the background mesh
  // and create two new edges at  boundary
  // if the first point of the shock  is at the boundary
  if (iedg1>0) { splitEdges();
                 createNewEdges(IP.at(0), ISH.at(0),IP.at(1), ISH.at(1)); }
}

//--------------------------------------------------------------------------//

void FxMshSps::setShockIndeces(unsigned nbDiscontinuities,unsigned ISPPNTS)
{
  unsigned I; 
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*r_SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*r_SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (r_nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//--------------------------------------------------------------------------//

void FxMshSps::setShockPointCoor(string zone, unsigned IP, unsigned ISH)
{
  if (zone=="Up")        { xsh = (*r_XYShu)(0,IP,ISH);
                           ysh = (*r_XYShu)(1,IP,ISH); }
  else if (zone=="Down") { xsh = (*r_XYShd)(0,IP,ISH);
                           ysh = (*r_XYShd)(1,IP,ISH); }
  logfile ("xsh: ", xsh, " ysh: ", ysh);
}

//--------------------------------------------------------------------------//

void FxMshSps::checkShockBndryEdgeCrossing()
{
  if(s1<0 || s1>1 || s2 <0 || s2>1) {
   logfile ("s(1): ", s1, " s(2): ", s2, "out of bonds\n");
   exit(1);
  }
  int iedge1 = iedg1+1; // only for file log printing, c++ indeces start from 0
  int iedge2 = iedg2+1; // only for file log printing, c++ indeces start from 0
  if (iedg2 != iedg1) {
   logfile ("Shock points (1) (2) NOT on the same boundary edge\n");
   logfile (iedge1, " ", iedge2, "\n");
   logfile (s1, " ", s2, "\n");
  }
}

//--------------------------------------------------------------------------//

void FxMshSps::splitEdges()
{
  I1 = (*bndfac)(0,iedg1); // c++ indeces start from 0
  I2 = (*bndfac)(1,iedg1); // c++ indeces start from 0 
  IBC = (*bndfac)(2,iedg1); // c++ indeces start from 0 

  (*bndfac)(2,iedg1) = -IBC;
  int iedg = iedg1+1; // only for file log printing, c++ indeces start from 0
  logfile("Removing background edge ", iedg, " ");
  logfile(I1, " ", I2, "\n");
}

//--------------------------------------------------------------------------//

void FxMshSps::createNewEdges(unsigned IP, unsigned ISH)
{ 
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistd(IP,ISH);
  }
  else {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistu(IP,ISH);
  }
  
  ibfac++;
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=ISHPlistu(IP,ISH);
   (*bndfac)(1,ibfac) = I2;
  }
  else {
   (*bndfac)(0,ibfac)=ISHPlistd(IP,ISH);
   (*bndfac)(1,ibfac) = I2;
  }
  ibfac++;
}

//--------------------------------------------------------------------------//

void FxMshSps::createNewEdges(unsigned IP1, unsigned ISH1,
                              unsigned IP2, unsigned ISH2)
{
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistd(IP2,ISH2);
  }
  else {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistu(IP1,ISH1);
  }

  ibfac++;
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=ISHPlistu(IP1,ISH1);
   (*bndfac)(1,ibfac) = I2;
  }
  else {
   (*bndfac)(0,ibfac)=ISHPlistd(IP2,ISH2);
   (*bndfac)(1,ibfac) = I2;
  }
  ibfac++;
}

//--------------------------------------------------------------------------//

void FxMshSps::setAddress()
{
  unsigned start;
  start = 0;
  r_XY = new Array2D <double> ((*ndim),(*npoin),&coor->at(start));
  start = (*npoin) * (*ndim);
  r_XYShu =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
  start = (*npoin) * (*ndim) + (*npshmax) * (*nshmax) * (*ndim);
  r_XYShd =
    new Array3D <double> ((*ndim),(*npshmax),(*nshmax),&coor->at(start));
}

//--------------------------------------------------------------------------//

void FxMshSps::setMeshData()
{
  nbfac = MeshData::getInstance().getData <unsigned> ("NBFAC");
  nbfacSh = MeshData::getInstance().getData <unsigned> ("NBFACSH");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem = MeshData::getInstance().getData <unsigned> ("NELEM");
  npoin = MeshData::getInstance().getData <unsigned> ("NPOIN");
  coor = MeshData::getInstance().getData <vector<double> > ("COOR");
  bndfac = MeshData::getInstance().getData <Array2D <int> > ("BNDFAC");
}

//--------------------------------------------------------------------------//

void FxMshSps::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  nshmax = PhysicsData::getInstance().getData <unsigned> ("NSHMAX");
  npshmax = PhysicsData::getInstance().getData <unsigned> ("NPSHMAX");
  r_nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  r_nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  r_nShockPoints = 
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  r_nShockEdges = 
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  r_typeSpecPoints = 
      PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  r_SHinSPPs = 
      PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
