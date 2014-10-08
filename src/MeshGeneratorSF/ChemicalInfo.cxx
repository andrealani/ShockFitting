// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ChemicalInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/IOFunctions.hh"
#include "Framework/MeshGenerator.hh"
#include "Framework/Log.hh"
#include "Common/PhysicsData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ChemicalInfo, MeshGenerator> readChemInfoProv("ChemicalInfo");

//--------------------------------------------------------------------------//

ChemicalInfo::ChemicalInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{
  m_Qref = 1;
  addOption("Qref",&m_Qref,
            "reference speed");
  m_IE = 1;
  addOption("IE",&m_IE,
             "Global index");
  m_IX = 1;
  addOption("IX",&m_IX,
             "Global index");
  m_IY = 1;
  addOption("IY",&m_IY,
             "Global index");
  m_IEV = 1;
  addOption("IEV",&m_IEV,
             "Global index");
  m_model = 1;
  addOption("MODEL",&m_model,
	    "Gas model");
  m_mixture = 1;
  addOption("MIXTURE",&m_mixture,
	    "Gas Mixture");
}

//--------------------------------------------------------------------------//

ChemicalInfo::~ChemicalInfo()
{
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setup()
{
  LogToScreen(VERBOSE, "ChemicalInfo::setup() => start\n");

  setPhysicsData();

  logfile.Open(getClassName());

  LogToScreen(VERBOSE, "ChemicalInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::unsetup()
{
  LogToScreen(VERBOSE, "ChemicalInfo::unsetup()\n");

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ChemicalInfo::generate()
{
  LogToScreen(INFO, "ChemicalInfo::generate()\n");

  setModel();
  setMixture();
  if (m_model == "PG") {
   (*nsp) = 1; (*ndof)=4;
   logfile("Number of variables: ",*ndof);
    return;}

  else if (m_model == "Cneq" || m_model == "TCneq") {
 
   file.open(getInputFiles().c_str());

   setMixtureFileName();
   setnbSpecies();
 
   setSize();
   setGlobalIndex();
   setnbDof();

   setSpeciesNames();
   setMolWeights();
   setFormEnthalp();
   setVibrTemp();
   setSpecHeatRatio();
   setMolTypes();

   file.close();
  }

  else {
   logfile("Model not implemented");
   exit(1);
  }
}

//--------------------------------------------------------------------------//

std::string ChemicalInfo::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setModel()
{
  model->at(0) = m_model;
  logfile("Thermodynamic model: ",m_model);
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setMixture() 
{
  mixture->at(0) = m_mixture;
  logfile("File mixture data: ",getInputFiles());
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setMixtureFileName()
{
  file >> dummystr;
  file >> mixtureFile_name;
  if (mixtureFile_name != m_mixture){
   logfile("error data file is wrong");
   logfile(mixtureFile_name, m_mixture.c_str());
   exit(1);}
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setnbSpecies()
{
  file >> dummystr;
  file >> *nsp;
  logfile( "Number of the species: ",*nsp);
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setSpeciesNames()
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> name->at(ISP);
   logfile("Chemical species: ",name->at(ISP));
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setMolWeights()
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> mm->at(ISP);
   logfile("MM: ",mm->at(ISP));
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setFormEnthalp() 
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> hf->at(ISP);
   logfile("HF: ",hf->at(ISP));
   if (mixture->at(0)=="ar4") { hf->at(ISP) = hf->at(ISP)/(m_Qref*m_Qref);}
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setVibrTemp()
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> thev->at(ISP);
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setSpecHeatRatio()
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> gams->at(ISP);
   logfile("GAMS: ",gams->at(ISP));
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setMolTypes()
{
  file >> dummystr;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> typemol->at(ISP);
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setSize()
{
  model->resize(1);
  mixture->resize(1);
  name->resize(*nsp);
  mm->resize(*nsp);
  hf->resize(*nsp);
  thev->resize(*nsp);
  gams->resize(*nsp);
  typemol->resize(*nsp);
} 

//--------------------------------------------------------------------------//

void ChemicalInfo::setGlobalIndex()
{
  *ie = m_IE + (*nsp); // = NSP+1
  *ix = m_IX + (*nsp); // = NSP+2
  *iy = m_IY + (*nsp); // = NSP+3
  *iev = m_IEV + (*nsp); // = NSP+4
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setnbDof()
{
  if (m_model=="Cneq") { (*ndof)=(*nsp)+(*ndim)+1; }
  else if (m_model =="TCneq") { (*ndof)=(*nsp)+(*ndim)+2; }
  else {
   logfile("Model not implemented");
   exit(1);
  }
  logfile("Number of variables: ",*ndof);
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setPhysicsData()
{
  ndim = PhysicsData::getInstance().getData <unsigned> ("NDIM");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  model = PhysicsData::getInstance().getData <std::vector<std::string> > ("MODEL");
  mixture = PhysicsData::getInstance().getData <std::vector<std::string> > ("MIXTURE");
  ie = PhysicsData::getInstance().getData <unsigned> ("IE");
  ix = PhysicsData::getInstance().getData <unsigned> ("IX");
  iy = PhysicsData::getInstance().getData <unsigned> ("IY");
  iev = PhysicsData::getInstance().getData <unsigned> ("IEV");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  name = PhysicsData::getInstance().getData <std::vector<std::string> > ("NAMESP");
  mm = PhysicsData::getInstance().getData <std::vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <std::vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <std::vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <std::vector<double> > ("GAMS");
  typemol = PhysicsData::getInstance().getData <std::vector<std::string> > ("TYPEMOL");
}
//--------------------------------------------------------------------------//

} //namespace ShockFitting
