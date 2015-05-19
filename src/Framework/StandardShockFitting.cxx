// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include <sstream>
#include "Framework/StandardShockFitting.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/MeshData.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/ConfigFileReader.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<StandardShockFitting, ShockFittingObj>
standardShockFittingProv("StandardShockFitting");

//--------------------------------------------------------------------------//

StandardShockFitting::StandardShockFitting(const std::string& objectName) :
  ShockFittingObj(objectName),
  m_createshockfile(),
  m_createTriangleFiles(),
  m_readInputFile1(),
  m_meshBackup(),
  m_readInputFile2(),
  m_bndryNodePtr(),
  m_bndryFacePtr(),
  m_redistrEqShockPoints(),
  m_findPhantPoints(),
  m_changeBndryPoints(),
  m_computeNormalVector(),
  m_computeShockLayer(),
  m_fixMeshSpecialPoints(),
  m_writeTriangleFile(),
  m_writeTriangleFileFreezedConnect(),
  m_callTriangle(),
  m_callTriangleLib(),
  m_triangleToCFfmt(),
  m_COOLFluiD(),
  m_CFmeshToTriangle(),
  m_copyZRoe1_0(),
  m_updateSolution(),
  m_fixSpecPoints(),
  m_copyZRoeSh0_1(),
  m_moveShPoints(),
  m_updatePhantPoints(),
  m_redistrShockPoints(),
  m_writeBackTriangleFile(),
  m_writeShockInfo(),
  m_meshRestore()
{
  m_version = "dummyVersion";
  addOption("Version",&m_version,
            "Current version of the StandardShockFitting");
  m_startFiles = false;
  addOption("startFromCapturedFiles",&m_startFiles,
            "The starting files are the ones from the captured solution");
  m_computeShockFittingResidual = false;
  addOption("shockFittingResidual",&m_computeShockFittingResidual,
            "Specifies if the shock fitting residual are computed");
}

//--------------------------------------------------------------------------//

StandardShockFitting::~StandardShockFitting()
{
}

//--------------------------------------------------------------------------//

void StandardShockFitting::configure(SConfig::OptionMap& cmap,
                                     const std::string& prefix)
{
  LogToScreen(VERBOSE, "StandardShockFitting::configure() => start\n");

  ShockFittingObj::configure(cmap, prefix);

  LogToScreen(VERBOSE, "StandardShockFitting::configure() => start\n");
}

//--------------------------------------------------------------------------//

void StandardShockFitting::setup()
{
  LogToScreen(VERBOSE, "StandardShockFitting::setup() => start\n");

  ShockFittingObj::setup();

  m_createshockfile = m_fConverter[0].ptr();
  m_createTriangleFiles = m_fConverter[1].ptr();
  m_readInputFile1 = m_mGenerator[0].ptr();
  m_meshBackup = m_cMaker[0].ptr();
  m_readInputFile2 = m_mGenerator[1].ptr();
  m_bndryNodePtr = m_fRemeshing[0].ptr();
  if(MeshData::getInstance().cellsFreezed()) {
   m_bndryFacePtr = m_fRemeshing[7].ptr(); }
  m_redistrEqShockPoints = m_fRemeshing[1].ptr();
  m_findPhantPoints = m_fRemeshing[2].ptr();
  m_changeBndryPoints = m_fRemeshing[3].ptr();
  m_computeNormalVector = m_cNormalVector.ptr();
  m_computeShockLayer = m_fRemeshing[4].ptr();
  m_fixMeshSpecialPoints = m_fRemeshing[5].ptr();
  m_writeTriangleFile = m_wMesh[0].ptr();
  if(MeshData::getInstance().freezedConnectivityOption()) {
   m_writeTriangleFileFreezedConnect = m_wMesh[3].ptr(); }
  m_callTriangle = m_mGenerator[2].ptr();
  m_callTriangleLib = m_mGenerator[3].ptr();
  if(m_computeShockFittingResidual) { 
   m_computeSFresidual = m_sUpdater[2].ptr();}
  m_triangleToCFfmt = m_fConverter[2].ptr();
  m_COOLFluiD = m_CFDSolver.ptr();
  m_CFmeshToTriangle = m_fConverter[3].ptr();
  m_copyZRoe1_0 = m_cMaker[1].ptr();
  m_updateSolution = m_cState.ptr();
  m_fixSpecPoints = m_sUpdater[0].ptr(); 
  m_copyZRoeSh0_1 = m_cMaker[2].ptr();
  m_moveShPoints = m_moveDps.ptr();
  m_updatePhantPoints = m_sUpdater[1].ptr();
  m_redistrShockPoints = m_fRemeshing[6].ptr();
  m_writeBackTriangleFile = m_wMesh[1].ptr();
  m_writeShockInfo = m_wMesh[2].ptr();
  m_meshRestore = m_cMaker[3].ptr();

  LogToScreen(VERBOSE, "StandardShockFitting::setup() => end\n");  
}

//--------------------------------------------------------------------------//

void StandardShockFitting::unsetup()
{
  LogToScreen(VERBOSE, "StandardShockFitting::unsetup() => start\n");

  ShockFittingObj::unsetup();
  
  LogToScreen(VERBOSE, "StandardShockFitting::unsetup() => end\n");
}

//--------------------------------------------------------------------------//

void StandardShockFitting::process()
{
  LogToScreen(VERBOSE, "StandardShockFitting::process() => start\n");

  // string for system command execution
  string execmd;

  // current file name
  stringstream* fname = MeshData::getInstance().getData <stringstream> ("FNAME");

  // name of the back file
  string* fnameback = MeshData::getInstance().getData <string> ("FNAMEBACK");

  ostringstream backdir;
  unsigned dummyIstep; 

  // set the Shock Fitting version
  MeshData::getInstance().setVersion(m_version);

  cout << "\n--------------------- Shock Fitting Solver ----------------------\n\n";
  cout << "_____________________ StandardShockFitting ______________________\n\n";

  cout << "StandardShockFitting.Version = " << m_version << "\n";
  cout << "_________________________________________________________________\n\n";

  cout << "            StandardShockFitting::pre-processing  \n";
  cout << "-----------------------------------------------------------------\n\n";

  cout << "Collecting general Physics informations \n\n";
 
  PhysicsData::getInstance().getPhysicsInfo()->read();
  PhysicsData::getInstance().getChemicalInfo()->read(); 
  PhysicsData::getInstance().getReferenceInfo()->read();
  cout << ".................................................\n";

  // the starting triangle files are generated from the captured solution.
  // A file format conversion is therefore askedd
  if(m_startFiles) {
   cout << "_________________________________________________\n\n";
   cout << "Creating starting SF files from the captured solution \n\n";

   m_createshockfile->convert();
   m_createTriangleFiles->convert();

   m_callTriangle->generate(string("na00.node"));
   cout << ".................................................\n";

   system(string("mv na00.poly na99.poly").c_str());
  }

  cout << "_________________________________________________\n\n";
  cout << "Building the initial computational domain\n\n";

  m_readInputFile1->generate();
  m_readInputFile2->generate();

  m_bndryNodePtr->remesh();

  if(MeshData::getInstance().cellsFreezed()) { m_bndryFacePtr->remesh(); }

  m_meshBackup->copy();

  m_redistrEqShockPoints->remesh();

  cout << ".................................................\n";

  cout << "_________________________________________________\n\n";


  cout << "\n-----------------------------------------------------------------";
  cout << "\n-----------------------------------------------------------------\n\n";
  cout << "              StandardShockFitting::starting the time loop   \n";
  cout << "-----------------------------------------------------------------\n";
  cout << "-----------------------------------------------------------------\n\n";

  for(unsigned I=MeshData::getInstance().getnbBegin();
    I<MeshData::getInstance().getnbSteps(); I++) {

   MeshData::getInstance().setIstep(I+1);

   cout << "              StandardShockFitting::step number => ";
   cout << MeshData::getInstance().getIstep() << "   \n";
   cout << "-----------------------------------------------------------------\n \n";

   if(MeshData::getInstance().getFreezedConnectivity()) {
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::";
    cout << "\n (!) StandardShockFitting::warning => freezed connectivity\n";
    cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n";
   }

   m_findPhantPoints->remesh();

   m_changeBndryPoints->remesh();

   m_computeNormalVector->remesh();

   m_computeShockLayer->remesh();
   m_fixMeshSpecialPoints->remesh();

   // if the connectivity is not freezed call triangle mesh generator
   if(!MeshData::getInstance().getFreezedConnectivity()) {
 
    if      (m_version=="original")  { m_writeTriangleFile->write();
                                      m_callTriangle->generate();    }
 
    else if (m_version=="optimized") { m_callTriangleLib->generate(); }
   }

   // if the connectivity is freezed do not call triangle mesh generator
   else if (MeshData::getInstance().getFreezedConnectivity()) {
    m_writeTriangleFileFreezedConnect->write();
   }

   m_triangleToCFfmt->convert();

   cout << "_________________________________________________________________\n\n";

   m_COOLFluiD->call();

   // change COOLFluiD output file name
   if(MeshData::getInstance().getnbProcessors()==1) {
    if(MeshData::getInstance().withP0()) { 
     execmd = "cp -f cfout-P0.CFmesh cfout.CFmesh"; 
     system(execmd.c_str());
     if(system(execmd.c_str())!=0) {
     cout << "StandardShockFitting::error => CFmesh file doesn't exist\n";
     exit(1); }
    }
    if(MeshData::getInstance().withP0()) { execmd = "rm -f cfout-P0.CFmesh"; 
                                           system(execmd.c_str()); }
   }
   else if (MeshData::getInstance().getnbProcessors()>1) {
    if(MeshData::getInstance().withP0()) { 
     execmd = "cp -f cfout-P?.CFmesh cfout.CFmesh";
     system(execmd.c_str());
     if(system(execmd.c_str())!=0) {
      cout << "StandardShockFitting::error => CFmesh file doesn't exist\n";
      exit(1); }
    }
    if(MeshData::getInstance().withP0()) { execmd = "rm -f cfout-P?.CFmesh";
                                           system(execmd.c_str()); }
   }

   cout << "_________________________________________________________________\n\n";

   m_CFmeshToTriangle->convert();

   if  (m_version=="original" )  { m_readInputFile1->generate(); }

   m_copyZRoe1_0->copy();

   m_updateSolution->update();
   m_fixSpecPoints->update();

   m_copyZRoeSh0_1->copy();

   m_moveShPoints->update();

   m_updatePhantPoints->update();

   if( I<1000 ) { m_redistrShockPoints->remesh(); }

   // if the solution must be saved in the I-step, the shock data
   // background grid data are written on output files
   if((I)%MeshData::getInstance().getnbIbak()==0) {
    m_writeBackTriangleFile->write();
    m_writeShockInfo->write();
   }

   m_meshRestore->copy();

   if(m_computeShockFittingResidual) { m_computeSFresidual->update(); }

   cout << "_________________________________________________________________\n";

   cout << "_________________________________________________________________\n\n";

   // create the directory to backup files
   backdir.str(string());
   unsigned nbDig=0;
   dummyIstep = I+1;
   while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
   backdir << setw(9-nbDig) << setfill('0') << left << string("step").c_str() << I+1;

   // during the current step the solution will be saved
   if((I)%MeshData::getInstance().getnbIbak()==0) {
    execmd = "mkdir " + backdir.str();
    system(execmd.c_str());

    execmd = "mv -f shocknor.dat sh99.dat cfout.CFmesh cfin.CFmesh "
             + *fnameback + ".node " + "cfin*plt "; 
    if (MeshData::getInstance().getVersion()=="original")
     { execmd = execmd + fname->str() + ".* ";}
    execmd = execmd + backdir.str();
    system(execmd.c_str());

    if (MeshData::getInstance().getnbProcessors()==1) {
     execmd = "mv cfout.plt cf" + backdir.str().substr(4,9) + ".plt";
     system(execmd.c_str());
     execmd = "mv cfout.surf.plt cf" + backdir.str().substr(4,9) + "-surf.plt";
     system(execmd.c_str());
     execmd = "mv -f cf*.plt " + backdir.str();
//     execmd = "cp Wall.plt-1 wall" + backdir.str().substr(4,9) + ".plt";
    }
    else if (MeshData::getInstance().getnbProcessors()>1) {
     execmd = "rename out cf"+backdir.str().substr(4,9)+" cfout-P?.plt";
     execmd = "mv -f cf*.plt " + backdir.str();
     system(execmd.c_str());
//     execmd = "cp Wall.plt-1 wall" + backdir.str().substr(4,9) + ".plt";
    }

    system(execmd.c_str());
   }

   // during the current step the solution wont be saved
   else {
    execmd = "rm -f shocknor.dat ";
    if(MeshData::getInstance().getVersion()=="original" && 
       !MeshData::getInstance().getFreezedConnectivity())
      { execmd = execmd + fname->str() +".* "; }
    system(execmd.c_str());
   }

   execmd = "cut -c1- residual.dat >> convergence.dat";
   system(execmd.c_str());
  }

  cout << "_________________________________________________________________\n";
  cout << "_________________________________________________________________\n";

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
