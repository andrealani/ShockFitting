// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/ComputeHeatFlux2DPg.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsConsts.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ComputeHeatFlux2DPg, StateUpdater> 
 computeHeatFlux2dPgProv("ComputeHeatFlux2DPg");

//--------------------------------------------------------------------------//

ComputeHeatFlux2DPg::ComputeHeatFlux2DPg(const std::string& objectName) :
  StateUpdater(objectName)
{

}

//--------------------------------------------------------------------------//

ComputeHeatFlux2DPg::~ComputeHeatFlux2DPg()
{
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::setup()
{
  LogToScreen(VERBOSE,"ComputeHeatFlux2DPg::setup() => start\n");
  
  m_paramToprimDimensional.ptr()->setup();
  
  LogToScreen(VERBOSE,"ComputeHeatFlux2DPg::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::unsetup()
{
  LogToScreen(VERBOSE,"ComputeHeatFlux2DPg::unsetup()\n");

  m_paramToprimDimensional.ptr()->unsetup();
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::configure(OptionMap& cmap, const std::string& prefix)
{
  StateUpdater::configure(cmap, prefix);

  // assign strings read in input.case to the object making the 
  // variable transformation
  m_paramToprimDimensional.name() = string("Param2PrimPg") + string("Dimensional");

  if (ConfigFileReader::isFirstConfig()) {
   m_paramToprimDimensional.ptr().reset
     (SConfig::Factory<VariableTransformer>::getInstance().
      getProvider(m_paramToprimDimensional.name())->
      create(m_paramToprimDimensional.name()));
  }

  // configure the object transforming the variables
  configureDeps(cmap, m_paramToprimDimensional.ptr().get());
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::update()
{
  LogToScreen(INFO,"ComputeHeatFlux2DPg::update()\n");

  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // working array
  vector<double> m_zroe(*ndof);
  vector<double> m_prim(*ndof);
  vector<double> m_XY(PhysicsInfo::getnbDim(),1);
  vector<double> distance(2);
  unsigned JVERT, KVERT, m_I;

  // vector storing temperature values
  vector<double> temperature(npoin->at(0));

  // array storing the IDnode of the cells that have one edge leading to the
  // "Wall" boundary
  Array2D<unsigned> wallcellnod(3,nbfac->at(0));

  // array storing heatflux components in the wall points
  // it will be resized after
  Array2D<double> heatflux(2,npoin->at(0));

  // vector storing the wall point of the wall cell in which the heat flux is computed
  // it will be resized after
  vector<unsigned> IDwallPoint(npoin->at(0));

  // working array storing gradient components in the wall points
  vector<double> gradT(2);

  // working variable storing the heatflux module
  double m_heatfluxMod;

  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {

   // assign the zroe of the IPOIN to the working vector m_zroe
   for(unsigned IV=0; IV<(*ndof); IV++) {
    m_zroe.at(IV) = (*zroe)(IV,IPOIN);
   }

   // convert the zroe variables in primitive variables
   // the XY coordinates are actually not used
   m_paramToprimDimensional.ptr()->transform(m_zroe,m_XY,m_prim);

   // assign the transformed prim variables of the IPOIN-point  to 
   // the corrisponding array
   temperature.at(IPOIN) = m_prim.at(3);
  }

  unsigned IWALLCELL=0;  
  for(unsigned IBFAC=0;IBFAC<nbfac->at(0);IBFAC++) {

   // this part must be changed using BCmap !!!
   if((*bndfac)(2,IBFAC) == 1) {

    wallcellnod(0,IWALLCELL) = (*bndfac)(0,IBFAC);
    wallcellnod(1,IWALLCELL) = (*bndfac)(1,IBFAC);

    // do a loop on the elements to look for the third node of the wall cell
    for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
     for(unsigned IVERT=0; IVERT<(*nvt); IVERT++) {
      if((*celnod)(IVERT,IELEM)==wallcellnod(0,IWALLCELL)) {
       if     (IVERT==0) {JVERT=1; KVERT=2;}
       else if(IVERT==1) {JVERT=0; KVERT=2;}
       else              {JVERT=0; KVERT=1;}
       if((*celnod)(JVERT,IELEM)==wallcellnod(1,IWALLCELL)) {
        wallcellnod(2,IWALLCELL)=(*celnod)(KVERT,IELEM);
       }
       else if((*celnod)(KVERT,IELEM)==wallcellnod(1,IWALLCELL)) {
        wallcellnod(2,IWALLCELL)=(*celnod)(JVERT,IELEM);
       }
      }
      else if((*celnod)(IVERT,IELEM)==wallcellnod(1,IWALLCELL)) {
       if     (IVERT==0) {JVERT=1; KVERT=2;}
       else if(IVERT==1) {JVERT=0; KVERT=2;}
       else              {JVERT=0; KVERT=1;}
       if((*celnod)(JVERT,IELEM)==wallcellnod(0,IWALLCELL)) {
        wallcellnod(2,IWALLCELL)=(*celnod)(KVERT,IELEM);
       }
       else if((*celnod)(KVERT,IELEM)==wallcellnod(0,IWALLCELL)) {
        wallcellnod(2,IWALLCELL)=(*celnod)(JVERT,IELEM);
       }
      }
     } // for nvt

    } // for nelem

    distance.assign(2,0);
    for(unsigned I=0;I<2;I++) {
     for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
      distance.at(I) += pow((*XY)(IDIM,wallcellnod(2,IWALLCELL)-1)-
                            (*XY)(IDIM,wallcellnod(I,IWALLCELL)-1),2);
     }
    }
    if(distance.at(0)<distance.at(1)) m_I = 0;
    else                              m_I = 1;

    IDwallPoint.at(IWALLCELL) = wallcellnod(m_I,IWALLCELL)-1;

    double deltaT = temperature.at(wallcellnod(2,IWALLCELL)-1)-
                    temperature.at(wallcellnod(m_I,IWALLCELL)-1);

    double deltaD = 0;
    for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
     deltaD += pow((*XY)(IDIM,wallcellnod(2,IWALLCELL)-1)-
                   (*XY)(IDIM,wallcellnod(m_I,IWALLCELL)-1),2);
    }

    deltaD = sqrt(deltaD);

    double angle = ((*XY)(1,wallcellnod(2,IWALLCELL)-1)-
                    (*XY)(1,wallcellnod(m_I,IWALLCELL)-1))/
                   ((*XY)(0,wallcellnod(2,IWALLCELL)-1)-
                    (*XY)(0,wallcellnod(m_I,IWALLCELL)-1));

    double XunitVector = cos(angle);
    double YunitVector = sin(angle);

    gradT.at(0) = deltaT/deltaD * XunitVector;
    gradT.at(1) = deltaT/deltaD * YunitVector;

    for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
     heatflux(IDIM,IDwallPoint.at(IWALLCELL)) =
      -getConductivity(temperature.at(wallcellnod(m_I,IWALLCELL)-1))*gradT.at(IDIM);
    }
    ++IWALLCELL;
   } // if wall edge
  } // for nbfac

  unsigned wallcells = IWALLCELL;

  IDwallPoint.resize(wallcells);

  // write the heat flux values on the output file
  FILE* writeHeatFlux;
  writeHeatFlux = fopen("heatflux.dat","w");

  fprintf(writeHeatFlux,"%s","TITLE = Heat Flux\n");
  fprintf(writeHeatFlux,"%s","VARIABLES = \"x0\" \"x1\" \"heatF\" \n");
  fprintf(writeHeatFlux,"%s","ZONE T = \"Wall heat flux\"\n");
  fprintf(writeHeatFlux,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
  fprintf(writeHeatFlux,"%s %u %s","I= ",wallcells,", J=1, K=1, ZONETYPE=Ordered\n");
  fprintf(writeHeatFlux,"%s","DATAPACKING=POINT\n");
  fprintf(writeHeatFlux,"%s","DT = (SINGLE, SINGLE, SINGLE)\n");
 
  for(unsigned IWALLCELL=0;IWALLCELL<wallcells;IWALLCELL++) {
   m_heatfluxMod = 0;
   for(unsigned IDIM=0;IDIM<PhysicsInfo::getnbDim();IDIM++) {
    fprintf(writeHeatFlux,"%28.14F",(*XY)(IDIM,IDwallPoint.at(IWALLCELL)));
    m_heatfluxMod += heatflux(IDIM,IDwallPoint.at(IWALLCELL))*
                     heatflux(IDIM,IDwallPoint.at(IWALLCELL));
   }
   m_heatfluxMod = sqrt(m_heatfluxMod);
   fprintf(writeHeatFlux,"%28.14E",m_heatfluxMod);
   fprintf(writeHeatFlux,"%s","\n");
  }

  fclose(writeHeatFlux);
  
  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

double ComputeHeatFlux2DPg::getDynViscosity(double absT)
{
  return (PhysicsConsts::SuthC1() * pow(absT,1.5) / 
          (absT + PhysicsConsts::SuthC2()));
}

//--------------------------------------------------------------------------//

double ComputeHeatFlux2DPg::getConductivity(double absT)
{
  return (getDynViscosity(absT)*PhysicsInfo::getCp()/PhysicsConsts::Pr());  
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::setAddress()
{
  unsigned start;
  zroe = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                             npoin->at(0), &zroeVect->at(0) );
  XY = new Array2D<double>(PhysicsInfo::getnbDim(),
                           npoin->at(0), &coorVect->at(0) );
  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  bndfac = new Array2D<int> (3,nbfac->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                                  PhysicsInfo::getnbShEdgesMax(),
                               &bndfacVect->at(0));
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::freeArray()
{
  delete zroe; delete XY;
  delete bndfac; delete celnod;
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem =  MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> > ("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> > ("CELNOD");
  BCmap = MeshData::getInstance().getData <vector<string> >("BoundariesMap");
}

//--------------------------------------------------------------------------//

void ComputeHeatFlux2DPg::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace shockFitting
