// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/ComputeShockPointsState.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Area.hh"
#include "MathTools/Jcycl.hh"
#include "MathTools/MinMax.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

ComputeShockPointsState::ComputeShockPointsState()
{
}

//--------------------------------------------------------------------------//

ComputeShockPointsState::~ComputeShockPointsState()
{
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::setup(double shLayerThick)
{
  cout << "     => ComputeShockPointsState::setup()\n";

  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // resize vectors and array
  setSize();
 
  // set the distance used to extract the upstream and downstream coordinates
  setShockLayerThickness(shLayerThick);

  logfile.Open(getClassName());
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::unsetup()
{
  cout << "     => ComputeShockPointsState::unsetup()\n";
  
  logfile.Close();

  // de-allocate dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::
 extractDownstreamAndUpstreamPoints(Array3D<double> normals)
{
  cout << "     => ComputeShockPointsState::extractDownstreamAndUpstreamPoints()\n";

  m_shockLayerThickness *= MeshData::getInstance().getDXCELL();

  vector<double> X(npoin->at(0));
  vector<double> Y(npoin->at(0));
  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
   X.at(IPOIN) = (*XY)(0,IPOIN);
   Y.at(IPOIN) = (*XY)(1,IPOIN);
  }

  const double minDomainXCoordinate = *min_element(X.begin(),X.end());
  const double maxDomainXCoordinate = *max_element(X.begin(),X.end());
  const double minDomainYCoordinate = *min_element(Y.begin(),Y.end());
  const double maxDomainYCoordinate = *max_element(Y.begin(),Y.end());

  double m_Xu, m_Yu;
  double m_Xd, m_Yd;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   unsigned m_shPoint = 0;

   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    m_Xu = (*XYSh)(0,ISHPOIN,ISH) +
           normals(0,ISHPOIN,ISH) * m_shockLayerThickness;
    m_Yu = (*XYSh)(1,ISHPOIN,ISH) +
           normals(1,ISHPOIN,ISH) * m_shockLayerThickness;
    m_Xd = (*XYSh)(0,ISHPOIN,ISH) -
           normals(0,ISHPOIN,ISH) * m_shockLayerThickness;
    m_Yd = (*XYSh)(1,ISHPOIN,ISH) -
           normals(1,ISHPOIN,ISH) * m_shockLayerThickness;

    // check if the extracted shock points are out of the domain
    // they are erased
    if(m_Xu>=minDomainXCoordinate && m_Xu<=maxDomainXCoordinate &&
       m_Yu>=minDomainYCoordinate && m_Yu<=maxDomainYCoordinate &&
       m_Xd>=minDomainXCoordinate && m_Xd<=maxDomainXCoordinate && 
       m_Yd>=minDomainYCoordinate && m_Yd<=maxDomainYCoordinate  )
    {
      XYShu(0,m_shPoint,ISH) = m_Xu; XYShu(1,m_shPoint,ISH) = m_Yu;
      XYShd(0,m_shPoint,ISH) = m_Xd; XYShd(1,m_shPoint,ISH) = m_Yd;
      (*XYSh)(0,m_shPoint,ISH) = (*XYSh)(0,ISHPOIN,ISH);
      (*XYSh)(1,m_shPoint,ISH) = (*XYSh)(1,ISHPOIN,ISH);
      ++m_shPoint;
    }
    else {
     cout << "     => ComputeShockPointsState:: (!) warning => IDpoint ";
     cout << ISHPOIN+1 << " is out of the domain\n";
     cout << "\n" <<  "         up= " << m_Xu << "," << m_Yu;
     cout << "; down= " << m_Xd << "," << m_Yd;
     cout << " has been erased\n\n";
    }
   }
   // update the number of shock points because some of them could have been erased
   nShockPoints->at(ISH) = m_shPoint;
   nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;
  }


  // plot downstream and upstream coordinates
  // it is useful to check if they are extracted within the shock layer
  FILE* plotState;
  stringstream fileName;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   fileName.str(string());
   fileName << "./log/UpstreamDownstreamCoor_ShNb_";
   fileName << ISH+1 << ".dat";
   plotState = fopen(fileName.str().c_str(),"w");

   fprintf(plotState,"%s","TITLE = Upstream and downstream coor\n");
   fprintf(plotState,"%s","VARIABLES = \"x0\" \"x1\"\n");
   fprintf(plotState,"%s","ZONE T = \"Up and Down coord\"\n");
   fprintf(plotState,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
   fprintf(plotState,"%s %u %s","I= ",nShockPoints->at(ISH),", J=2, K=1, ZONETYPE=Ordered\n");
   fprintf(plotState,"%s","DATAPACKING=POINT\n");
   fprintf(plotState,"%s","DT = (SINGLE, SINGLE)\n");

   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    fprintf(plotState,"%22.14F",XYShu(0,ISHPOIN,ISH));
    fprintf(plotState,"%22.14F",XYShu(1,ISHPOIN,ISH));
    fprintf(plotState,"%s","\n");

   }

   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
   fprintf(plotState,"%22.14F",XYShd(0,ISHPOIN,ISH));
    fprintf(plotState,"%22.14F",XYShd(1,ISHPOIN,ISH));
    fprintf(plotState,"%s","\n");
   }

   fclose(plotState);
  }
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::interpDownstreamAndUpstreamState()
{
   cout << "     => ComputeShockPointsState::interpDownstreamAndUpstreamState()\n"; 

  // working variables and vectors
  double s, t, help;
  unsigned I;
  vector<double> xp(4);
  vector<double> yp(4);
  vector<unsigned> idxs(3);
  vector<double> a(3);

  // create object computing area
  Area boundedArea;

  // create object Jcycl
  Jcycl J;

  // create object finding minimum and maximum values
  MinMax <double> m;

  // boolean variables evaluating if the shock point has been located
  // in the mesh connectivity
  bool found;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {

    found = false;
    for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
     for(unsigned IV=0; IV<(*nvt); IV++) {
      I = (*celnod)(IV,IELEM); // global code number
      xp.at(IV) = (*XY)(0,I-1);
      yp.at(IV) = (*XY)(1,I-1);
     }
     xp.at(3) = XYShd(0,ISHPOIN,ISH);  // node to be located
     yp.at(3) = XYShd(1,ISHPOIN,ISH);  // node to be located

     idxs.at(0) = 0; idxs.at(1) = 1; idxs.at(2) = 2; // c++ indeces start from 0
     boundedArea.computeArea(xp,yp,idxs);

     help = 1/boundedArea.getArea();
     idxs.at(2) = 3; // node to be located

     // compute area coordinates (in the x-y plane)
     for(unsigned IV=0; IV<(*nvt); IV++) {
      idxs.at(0) = J.callJcycl(1+IV+1)-1;
      idxs.at(1) = J.callJcycl(2+IV+1)-1;
      boundedArea.computeArea(xp,yp,idxs);
      a.at(IV) = boundedArea.getArea() * help;
     }

     s = m.min(a.at(0), a.at(1), a.at(2));
     t = m.max(a.at(0), a.at(1), a.at(2));

     if (( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)) {

      for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) { (*ZRoeShd)(IDOF,ISHPOIN,ISH) = 0;}
      for(unsigned IV=0; IV<(*nvt); IV++) {
       I = (*celnod)(IV,IELEM);
       for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
        (*ZRoeShd)(IDOF,ISHPOIN,ISH) =
          (*ZRoeShd)(IDOF,ISHPOIN,ISH) + a.at(IV) * (*zroe)(IDOF,I-1);
       }
      }
      found = true; break;
     } // if ( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)
    } // for IELEM<nelem->at(0)

    if(found==false) {
     cout << "        ComputeShockPointsState::error => Search failed for shock point:\n";
     cout << "                                  " << XYShd(0,ISHPOIN,ISH);
     cout << " " << XYShd(1,ISHPOIN,ISH) << endl;
     cout << "\n  -> Try to reduce the shock layer thickness\n";
     exit(1);
    } // if found=false

    found = false;
    for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
     for(unsigned IV=0; IV<(*nvt); IV++) {
      I = (*celnod)(IV,IELEM); // global code number
      xp.at(IV) = (*XY)(0,I-1);
      yp.at(IV) = (*XY)(1,I-1);
     }
     xp.at(3) = XYShu(0,ISHPOIN,ISH);  // node to be located
     yp.at(3) = XYShu(1,ISHPOIN,ISH);  // node to be located

     idxs.at(0) = 0; idxs.at(1) = 1; idxs.at(2) = 2; // c++ indeces start from 0
     boundedArea.computeArea(xp,yp,idxs);
     help = 1/boundedArea.getArea();
     idxs.at(2) = 3; // node to be located

     // compute area coordinates (in the x-y plane)
     for(unsigned IV=0; IV<(*nvt); IV++) {
      idxs.at(0) = J.callJcycl(1+IV+1)-1;
      idxs.at(1) = J.callJcycl(2+IV+1)-1;
      boundedArea.computeArea(xp,yp,idxs);
      a.at(IV) = boundedArea.getArea() * help;
     }

     s = m.min(a.at(0), a.at(1), a.at(2));
     t = m.max(a.at(0), a.at(1), a.at(2));

     if (( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)) {

      for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) { (*ZRoeShu)(IDOF,ISHPOIN,ISH) = 0;}
      for(unsigned IV=0; IV<(*nvt); IV++) {
       I = (*celnod)(IV,IELEM);
       for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
        (*ZRoeShu)(IDOF,ISHPOIN,ISH) =
          (*ZRoeShu)(IDOF,ISHPOIN,ISH) + a.at(IV) * (*zroe)(IDOF,I-1);
       }
      }
      found = true; break;
     } // if ( s>= 0 && s<= 1) && ( t>= 0 && t<= 1)

    } // for IELEM<nelem->at(0)
    if(found==false) {
     cout << "        ComputeShockPointsState::error => Search failed for shock points:\n";
     cout << "                                  " << XYShu(0,ISHPOIN,ISH);
     cout << " " << XYShu(1,ISHPOIN,ISH) << endl;
     cout << "\n  -> Try to reduce the shock layer thickness\n";
     exit(1);
    } // if found=false

   } // for ISHPOIN<nShockPoints->at(ISH)
  } // for ISH<nShocks

/*  for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
(*ZRoeShu)(IDOF,0,0) = (*ZRoeShu)(IDOF,1,0);
(*ZRoeShd)(IDOF,0,0) = (*ZRoeShd)(IDOF,1,0);
(*ZRoeShd)(IDOF,nShockPoints->at(0)-1,0) = (*ZRoeShd)(IDOF,nShockPoints->at(0)-2,0);
(*ZRoeShu)(IDOF,nShockPoints->at(0)-1,0) = (*ZRoeShu)(IDOF,nShockPoints->at(0)-2,0);
  }
*/
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::assignDownstreamAndUpstreamState()
{
  cout << "     => ComputeShockPointsState::assignDownstreamAndUpstreamState()\n";

  logfile("Assign downstream and upstream state\n");

  logfile("Number of shocks: ",(*nShocks),"\n");

  // define varID variable used to point the parameter vector variables
  // that assign upstream and downstream state 
  // varID is chosen according to the gas model 
  unsigned varID;
  if(ChemicalInfo::getModel()=="PG") { varID = 1; }
  else if (ChemicalInfo::getModel()=="TCneq") { varID = (*nsp)+1; }
  else {
   cout << "        ComputeShockPointsState::error => gas model not implemented\n";
   exit(1);
  }

  double dum;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
   logfile("\n------------------------\n");
   logfile("Shock nb => ",ISH+1,"\n");
   logfile("number of shock points is => ",nShockPoints->at(ISH),"\n\n");

   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {

    logfile("Shock point: ",ISHPOIN+1);
    logfile(" -> upstream = ",(*ZRoeShu)(varID,ISHPOIN,ISH));
    logfile(", downstream = ",(*ZRoeShd)(varID,ISHPOIN,ISH),"\n");
    if( (*ZRoeShu)(varID,ISHPOIN,ISH)>(*ZRoeShd)(varID,ISHPOIN,ISH)) {
     logfile("(!) downstream and upstream state are inverted\n");
     for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
      dum = (*ZRoeShu)(IDOF,ISHPOIN,ISH);
      (*ZRoeShu)(IDOF,ISHPOIN,ISH) = (*ZRoeShd)(IDOF,ISHPOIN,ISH);
      (*ZRoeShd)(IDOF,ISHPOIN,ISH) = dum;
     }
    }

    logfile("Upstream:: ");
    for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
     logfile((*ZRoeShu)(IDOF,ISHPOIN,ISH), " ");
     // make a backup of the zroesh
     (*ZRoeShuOld)(IDOF,ISHPOIN,ISH) = (*ZRoeShu)(IDOF,ISHPOIN,ISH);
    }

    logfile("\nDownstream:: ");
    for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
     // make a backup of the zroesh
     logfile((*ZRoeShd)(IDOF,ISHPOIN,ISH), " ");
     (*ZRoeShdOld)(IDOF,ISHPOIN,ISH) = (*ZRoeShd)(IDOF,ISHPOIN,ISH);
    }
    logfile("\n");
   }
  }
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::setSize()
{
  XYShu.resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               (*nShocks));
  XYShd.resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               (*nShocks));
  ZRoeShuOld->resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  ZRoeShdOld->resize((*ndof),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XY = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0));
  zroe = new Array2D<double> (PhysicsInfo::getnbDofMax(),
                              totsize,
                              &zroeVect->at(0));
  celnod = new Array2D<int> ((*nvt),nelem->at(0),
                             &celnodVect->at(0));

  unsigned start = npoin->at(0)*PhysicsInfo::getnbDofMax();
  ZRoeShu =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::freeArray()
{
  delete XY; delete zroe; delete celnod;
  delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::setMeshData()
{
  nvt = MeshData::getInstance().getData<unsigned>("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void ComputeShockPointsState::setPhysicsData()
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
  ZRoeShuOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZRoeShdOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
