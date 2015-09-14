// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/DetectorAlgorithm.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "MathTools/LeastSquaresMethod.hh"
#include "SConfig/ObjectProvider.hh"
#include "ShockDetectorSF/ComputeShockPointsState.hh"
#include "ShockDetectorSF/FindShockSurfaceNormals2D.hh"
#include "ShockDetectorSF/IdentifySpecialPoints.hh"
#include "ShockDetectorSF/RedistributeShockPoints.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------/

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DetectorAlgorithm, ShockDetector>
detectorAlgorithmProv("DetectorAlgorithm");

//--------------------------------------------------------------------------//

DetectorAlgorithm::DetectorAlgorithm(const std::string& objectName) :
  ShockDetector(objectName)
{
  // you can choose between
  // 1) GnoffoShockSensor
  // 2) NormalMachNumberPg2D
  // 3) DensityGradientMaximaPg2D
  // GnoffoshockSensor works very well for the circular cylinder
  // the other two do not work very well (maybe it is because of the
  // evaluation of the gradient)
  m_whichDetector.name() = "";
  addOption("Detector",&m_whichDetector,
            "Specifies the detector algorithm");
  m_param2prim.name() = "dummyVariableTransformer";
  // you can choose between "Ellipse", "SplittingCurves" and "Polynomial"
  m_fittingTechnique = "";
  addOption("fittingTechnique",&m_fittingTechnique,
             "Set the technique used to fit the evaluated shock points");
  // if "Polynomial" is chosen the order must be specified
  m_polynomialOrder = 0;
  addOption("polynomialOrder",&m_polynomialOrder,
             "Set the polynomial order");
  // if "SplittingCurves" is chosen the number of x and y segments in which the
  // curve will be divided must be specified in addition to the polynomial order
  // to apply for each part of the splitted curve
  m_nbSegments = vector<unsigned>();
  addOption("nbXandYsegments",&m_nbSegments,
             "Specify the number of x and y segments in which the curve will be divided into");
  m_iSegPolynomialOrders = vector<unsigned>();
  addOption("segmPolynomialOrders",&m_iSegPolynomialOrders,
             "Specify the polynomial orders for each part of the splitted curve");
  // if "SplittingCurves" is chosen the smoothing option must be specified
  // default is set to false
  m_smoothOption = false;
  addOption("smoothingOption",&m_smoothOption,
             "Specify if the Y-ccordinates are smoothed after the leats squares");
  // set the shock layer thickness that is the distance used to extract the
  // upstream and downstream variables. 
  // It will be multiplied for the mesh characteristical length
  m_shockLayerThickness = 1.;
  addOption("shockLayerFactor",&m_shockLayerThickness,
             "Set the distance to extract the upstream and downstream coordinates"); 
}

//--------------------------------------------------------------------------/

DetectorAlgorithm::~DetectorAlgorithm()
{
}

//--------------------------------------------------------------------------//

void DetectorAlgorithm::setup()
{
  LogToScreen (VERBOSE, "DetectorAlgorithm::setup => start\n");

  m_whichDetector.ptr()->setup();

  m_param2prim.ptr()->setup();

  LogToScreen (VERBOSE, "DetectorAlgorithm::setup => end\n");
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::unsetup()
{
  LogToScreen (VERBOSE, "DetectorAlgorithm::unsetup()\n");
  
  m_whichDetector.ptr()->unsetup();

  m_param2prim.ptr()->unsetup();
}

//----------------------------------------------------------------------------//


void DetectorAlgorithm::configure(OptionMap& cmap, const std::string& prefix)
{
  ShockDetector::configure(cmap, prefix);

  // assign strings read on input.case file to variable transformer object
  m_param2prim.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;

  if (ConfigFileReader::isFirstConfig()) {
   m_param2prim.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                             getProvider(m_param2prim.name())
                             ->create(m_param2prim.name()));

   m_whichDetector.ptr().reset(SConfig::Factory<ShockDetector>::getInstance().getProvider(m_whichDetector.name())->create(m_whichDetector.name()));
  }

  // configure object detecting discontinuities
  configureDeps(cmap, m_whichDetector.ptr().get());

  // configure object transforming values
  configureDeps(cmap, m_param2prim.ptr().get());
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::detect()
{
  LogToScreen(INFO,"DetectAlgorithm::using => "+
                       m_whichDetector.name()+"\n");

  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // create object fitting the shock points distribution
  LeastSquaresMethod LSM;

  // create object re-distribute shock points
  RedistributeShockPoints newShPointsDistribution;

  // create object locating special points
  IdentifySpecialPoints mySpecialPoints;

  // create object compute normal vector to the shock points
  FindShockSurfaceNormals2D shockSurfNormals;

  // create object evaluate upstream and downstream state for each
  // shock point
  ComputeShockPointsState shockPointsUpAndDownState;

  // working array
  vector<double> m_zroe(*ndof);
  vector<double> m_prim(*ndof);
  vector<double> m_XY(PhysicsInfo::getnbDim(),1);

  // plot the backgournd grid
  FILE* gridfile;
  gridfile = fopen("log/BackgroundGrid.dat","w");

  stringstream var;
  var.str(string());
  for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
   var << "\"Z[" << IDOF << "]\" ";
  }

  fprintf(gridfile,"%s","TITLE = Background Grid\n");
  fprintf(gridfile,"%s %s %s","VARIABLES = \"x0\" \"x1\"",var.str().c_str(), "\n");
  fprintf(gridfile,"%s %u","ZONE N=",npoin->at(0));
  fprintf(gridfile,"%s %u"," T= \"Grid\", E = ", nelem->at(0));
  fprintf(gridfile,"%s",", F=FEPOINT, ET=TRIANGLE, SOLUTIONTIME=0\n");

  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
   for(unsigned IV=0;IV<PhysicsInfo::getnbDim();IV++) {
    fprintf(gridfile,"%20.14F %s",(*XY)(IV,IPOIN)," ");
   }
   for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
    fprintf(gridfile,"%20.14F %s",(*zroe)(IDOF,IPOIN)," ");
   }
  fprintf(gridfile,"%s","\n");
  }

  for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
   for(unsigned IVERT=0;IVERT<(*nvt);IVERT++) {
    fprintf(gridfile,"%11u",(*celnod)(IVERT,IELEM));
   }
   fprintf(gridfile,"%s","\n");
  }
  fclose(gridfile);

  // convert parameter vector variables in primitive variables
  // in order to use the shock sensor
  prim.resize(npoin->at(0)*(*ndof));
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {

   // assign the zroe of the IPOIN to the working vector m_zroe
   for(unsigned IV=0; IV<(*ndof); IV++) {
    m_zroe.at(IV) = (*zroe)(IV,IPOIN);
   }

   // convert the zroe variables in primitive variables
   // the XY coordinates are actually not used
   m_param2prim.ptr()->transform(m_zroe,m_XY,m_prim);
 
   // assign the transformed prim variables of the IPOIN-point  to 
   // the corrisponding array
   for(unsigned IV=0; IV<(*ndof); IV++) {
    prim.at(IPOIN*(*ndof)+IV) = m_prim.at(IV);
   }
  }
           
  m_whichDetector.ptr()->detect(prim);

  // the shocks points coordinates are stored as if there is 
  // only one shock. The number of shocks will be evaluated after
  unsigned ISH=0;

  if(m_fittingTechnique!="") {

   // define working vectors storing shock points coordinates
   std::vector<double> XSh(nShockPoints->at(ISH));
   std::vector<double> YSh(nShockPoints->at(ISH));

   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    XSh.at(ISHPOIN) = (*XYSh)(0,ISHPOIN,ISH);
    YSh.at(ISHPOIN) = (*XYSh)(1,ISHPOIN,ISH);
   }
   if(m_fittingTechnique=="Ellipse") 
    { LSM.fitEllipse(XSh,YSh); }
   else if (m_fittingTechnique=="SplittingCurves") 
    {LSM.fitSplittingCurves(XSh,YSh,
                            m_nbSegments,m_iSegPolynomialOrders, m_smoothOption); }
   else if (m_fittingTechnique=="Polynomial") {
    if(m_polynomialOrder==0) {
     cout << "DetectorAlgorithm::error => the polynomial order is not specified\n";
     exit(1);
    }
    else { LSM.fitData(XSh,YSh,m_polynomialOrder); }
   }
   else {
    cout << "DetectorAlgorithm::error => ";
    cout << m_fittingTechnique << " fitting technique not implemented\n";
    exit(1);
   }

   XSh.resize(LSM.getNewPoints());
   YSh.resize(LSM.getNewPoints());

   XSh = LSM.getXvector();
   YSh = LSM.getYvector();

   nShockPoints->at(ISH) = LSM.getNewPoints();
   nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;

   // sort the shock points starting from the one in the top of the domain
   // and ending with the one on the bottom
 
   // first find the shock point on the top of the domain
   unsigned IDtopShPoint=0; 
   double maxYSh = YSh.at(IDtopShPoint);
   for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    if(YSh.at(ISHPOIN)>maxYSh) { IDtopShPoint=ISHPOIN;
                                 maxYSh = YSh.at(IDtopShPoint);}
   } 

   // compute the distance between the shock points and the
   // shock point on the top of the domain
   vector<unsigned> newIdISHPOIN(nShockPoints->at(ISH));
   vector<double> S(nShockPoints->at(ISH));
   for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
    for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
     S.at(ISHPOIN)= sqrt(pow(XSh.at(ISHPOIN)-XSh.at(IDtopShPoint),2)+
                         pow(YSh.at(ISHPOIN)-YSh.at(IDtopShPoint),2));
     // initialize the new ID used to sort the shock points
     newIdISHPOIN.at(ISHPOIN) = ISHPOIN;
    }
   }

   // sort the new ID shock points vector using the distance from the point
   // on the top of the domain
   sort(newIdISHPOIN.begin(),newIdISHPOIN.end(),
        ShockFitting::DetectorAlgorithm::indexSort(S));

   // initialize NODCODSH which is part of NODCOD
   // If the code -99 is used this means no shock point
   // If the code 10 is used this means shock point
   for (unsigned ISH=0; ISH < PhysicsInfo::getnbShMax(); ISH++) {
    for (unsigned K=0; K < PhysicsInfo::getnbShPointsMax(); K++)
     {(*nodcodSh)(K,ISH) = -99;}
  }

  // assign the new shock points coordinates using the vector
  // storing the new ID for the sorted shock points
  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
    for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
     (*XYSh)(0,ISHPOIN,ISH) = XSh.at(newIdISHPOIN.at(ISHPOIN));
     (*XYSh)(1,ISHPOIN,ISH) = YSh.at(newIdISHPOIN.at(ISHPOIN));
     (*nodcodSh)(ISHPOIN,ISH) = 10;    
    }
   }
  }

  mySpecialPoints.identify();

  // re-distribute the shock points
  newShPointsDistribution.distribute();
 
  // plot shock polynline
  FILE* shockPolyLine;
  shockPolyLine = fopen("log/ShockPolyline.dat","w");

  fprintf(shockPolyLine,"%s","TITLE = Shock polyline\n");
  fprintf(shockPolyLine,"%s","VARIABLES = \"x0\" \"x1\"\n");
  fprintf(shockPolyLine,"%s","ZONE T = \"ShockPolyline\"\n");
  fprintf(shockPolyLine,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
  fprintf(shockPolyLine,"%s %u %s","I= ",nShockPoints->at(ISH),", J=1, K=1, ZONETYPE=Ordered\n");
  fprintf(shockPolyLine,"%s","DATAPACKING=POINT\n");
  fprintf(shockPolyLine,"%s","DT = (SINGLE, SINGLE)\n");

  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
   fprintf(shockPolyLine,"%22.14F",(*XYSh)(0,ISHPOIN,ISH));
   fprintf(shockPolyLine,"%22.14F",(*XYSh)(1,ISHPOIN,ISH));
   fprintf(shockPolyLine,"%s","\n");
  }

  fclose(shockPolyLine);

  // compute the unit normal to the shock points
  shockSurfNormals.computeUnitNormals();

  // setup the vectors and array used in shockPointsUpAndDownState object
  shockPointsUpAndDownState.setup(m_shockLayerThickness);

  // extract the upstream and downstream points along the
  // normal vector of each shock point
  shockPointsUpAndDownState.extractDownstreamAndUpstreamPoints
   (shockSurfNormals.getUnitNormalsVector());

  // interpolate the downstream and upstream state for each shock point
  // using the state of the surrounding nodes
  shockPointsUpAndDownState.interpDownstreamAndUpstreamState();

  // assign the downstream and upstream state according to the value
  // of the interpolated variables
  shockPointsUpAndDownState.assignDownstreamAndUpstreamState();

  // de-allocate the vectors and array used in shockPointsUpAndDownState 
  // object
  shockPointsUpAndDownState.unsetup();

  // de-allocate dynamic arrays
  freeArray();
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::detect(vector<double>& firstInfoVector)
{
  LogToScreen(INFO,"DetectAlgorithm::detect()\n");
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::setAddress()
{
  unsigned totsize0 = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  zroe = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                               totsize0, &zroeVect->at(0) );
  XY = new Array2D<double>(PhysicsInfo::getnbDim(),
                           totsize0, &coorVect->at(0) );
  celnod = new Array2D<int>((*nvt),nelem->at(0),&celnodVect->at(0));
  nodcodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(npoin->at(0)));
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::freeArray()
{
  delete XY; delete zroe;
  delete celnod; delete nodcodSh;
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM"); 
  nodcod = MeshData::getInstance().getData <vector<int> > ("NODCOD");
  zroeVect =  MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//----------------------------------------------------------------------------//

void DetectorAlgorithm::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
   PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  XYSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
