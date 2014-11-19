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
  m_readInputFile1(),
  m_meshBackup(),
  m_readInputFile2(),
  m_bndryNodePtr(),
  m_redistrEqShockPoints(),
  m_findPhantPoints(),
  m_changeBndryPoints(),
  m_computeNormalVector(),
  m_computeShockLayer(),
  m_fixMeshSpecialPoints(),
  m_writeTriangleFile(),
  m_callTriangle(),
  m_triangleToCFmesh(),
  m_COOLFluiD(),
  m_CFmeshToTriangle(),
  m_readNewMesh(),
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


  m_readInputFile1 = m_mGenerator[0].ptr();
  m_meshBackup = m_cMaker[0].ptr();
  m_readInputFile2 = m_mGenerator[1].ptr();
  m_bndryNodePtr = m_fRemeshing[0].ptr();
  m_redistrEqShockPoints = m_fRemeshing[1].ptr();
  m_findPhantPoints = m_fRemeshing[2].ptr();
  m_changeBndryPoints = m_fRemeshing[3].ptr();
  m_computeNormalVector = m_cNormalVector.ptr();
  m_computeShockLayer = m_fRemeshing[4].ptr();
  m_fixMeshSpecialPoints = m_fRemeshing[5].ptr();
  m_writeTriangleFile = m_wMesh[0].ptr();
  m_callTriangle = m_mGenerator[2].ptr();
  m_triangleToCFmesh = m_fConverter[0].ptr();
  m_COOLFluiD = m_CFDSolver.ptr();
  m_CFmeshToTriangle = m_fConverter[1].ptr();
  m_readNewMesh = m_mGenerator[3].ptr();
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

  cout << "\n---------------- Shock Fitting Solver ----------------\n\n";
  cout << "________________ StandardShockFitting ________________\n\n";

  PhysicsData::getInstance().getPhysicsInfo()->read();
  PhysicsData::getInstance().getChemicalInfo()->read(); 
  PhysicsData::getInstance().getReferenceInfo()->read();

  m_readInputFile1->generate();

  m_readInputFile2->generate();

  m_bndryNodePtr->remesh();

  m_meshBackup->copy();

  m_redistrEqShockPoints->remesh();

  cout << "\n______________________________________________________\n\n";
  cout << "StandardShockFitting::entering in the first step";
  cout << "\n______________________________________________________\n\n";

  for(unsigned I=MeshData::getInstance().getnbBegin();
    I<MeshData::getInstance().getnbSteps(); I++) {

   MeshData::getInstance().setIstep(I+1);

   cout << "StandardShockFitting::step number => ";
   cout << MeshData::getInstance().getIstep() << endl;
   cout << "______________________________________________________\n\n";

   m_findPhantPoints->remesh();
   m_changeBndryPoints->remesh();
   m_computeNormalVector->remesh();
   m_computeShockLayer->remesh();
   m_fixMeshSpecialPoints->remesh();

   m_writeTriangleFile->write();

   m_callTriangle->generate();

   m_triangleToCFmesh->convert();

   cout << "______________________________________________________\n\n";

   m_COOLFluiD->call();

   // change COOLFluiD output file name
   if(MeshData::getInstance().getnbProcessors()==1) {
    execmd = "cp -f cfout-P0.CFmesh cfout.CFmesh"; system(execmd.c_str());
    if(system(execmd.c_str())!=0) {
    cout << "StandardShockFitting::error => CFmesh file doesn't exist\n";
    exit(1); }
    execmd = "rm -f cfout-P0.CFmesh"; system(execmd.c_str());
   }
   else if (MeshData::getInstance().getnbProcessors()>1) {
    execmd = "cp -f cfout-P?.CFmesh cfout.CFmesh"; system(execmd.c_str());
    if(system(execmd.c_str())!=0) {
     cout << "StandardShockFitting::error => CFmesh file doesn't exist\n";
     exit(1); }
    execmd = "rm -f cfout-P?.CFmesh"; system(execmd.c_str());
   }

   cout << "______________________________________________________\n\n";

   m_CFmeshToTriangle->convert();

   m_readNewMesh->generate();

   m_copyZRoe1_0->copy();

   m_updateSolution->update();

   m_fixSpecPoints->update();

   m_copyZRoeSh0_1->copy();

   m_moveShPoints->update();

   m_updatePhantPoints->update();

   m_redistrShockPoints->remesh();

   m_writeBackTriangleFile->write();
   m_writeShockInfo->write();

   m_meshRestore->copy();

   cout << "______________________________________________________\n\n";

   // create the directory to backup files
   backdir.str(string());
   unsigned nbDig=0;
   dummyIstep = I+1;
   while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
   backdir << setw(9-nbDig) << setfill('0') << left << string("step").c_str() << I+1;

   if((I)%MeshData::getInstance().getnbIbak()==0) {
    execmd = "mkdir " + backdir.str();
    system(execmd.c_str());

    execmd = "mv -f shocknor.dat sh99.dat cfout.CFmesh " + fname->str() + ".* "
             + *fnameback + ".node " + backdir.str();
    system(execmd.c_str());

    if (MeshData::getInstance().getnbProcessors()==1) {
     execmd = "cp -f cfout-P0.plt cf" + backdir.str().substr(4,9) + ".plt";
    }
    else if (MeshData::getInstance().getnbProcessors()>1) {
     execmd = "rename out cf"+backdir.str().substr(4,9)+" cfout-P?.plt";
     execmd = "mv -f cf*.plt " + backdir.str();
    }

    system(execmd.c_str());
   }

   else {
    execmd = "rm -f shocknor.dat "+ fname->str() +".* ";
    execmd = execmd + *fnameback+".node sh99.dat ";
    system(execmd.c_str());
   }

   execmd = "cut -c1- residual.dat >> convergenza.dat";
   system(execmd.c_str());

  }

  cout << "______________________________________________________\n";
  cout << "______________________________________________________\n";

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
