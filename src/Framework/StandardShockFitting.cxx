// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

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

  PhysicsData::getInstance().getPhysicsInfo()->read();
  PhysicsData::getInstance().getChemicalInfo()->read(); 
  PhysicsData::getInstance().getReferenceInfo()->read();

  m_readInputFile1->generate();

  m_readInputFile2->generate();

  m_bndryNodePtr->remesh();

  m_meshBackup->copy();

  m_redistrEqShockPoints->remesh();
  m_findPhantPoints->remesh();
  m_changeBndryPoints->remesh();
  m_computeNormalVector->remesh();
  m_computeShockLayer->remesh();
  m_fixMeshSpecialPoints->remesh();

  m_writeTriangleFile->write();

  m_callTriangle->generate();

  m_triangleToCFmesh->convert();

  m_COOLFluiD->call();

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

  LogToScreen(VERBOSE, "StandardShockFitting::process() => end\n");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
