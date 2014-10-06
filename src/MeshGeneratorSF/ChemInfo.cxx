// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ChemInfo.hh"
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
ObjectProvider<ChemInfo, MeshGenerator> readChemInfoProv("ChemInfo");

//--------------------------------------------------------------------------//

ChemInfo::ChemInfo(const std::string& objectName) :
  MeshGenerator(objectName)
{
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

ChemInfo::~ChemInfo()
{
}

//--------------------------------------------------------------------------//

void ChemInfo::setup()
{
  LogToScreen(VERBOSE, "ChemInfo::setup() => start\n");

  setPhysicsData();

  LogToScreen(VERBOSE, "ChemInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ChemInfo::unsetup()
{
  LogToScreen(VERBOSE, "ChemInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ChemInfo::generate()
{
  LogToScreen(INFO, "ChemInfo::generate()\n");

  setModel();
  setMixture();
 std::cout << "set model" << std::endl;
  if (m_model == "PG") {
   (*nsp) = 1; (*ndof)=4;
   //re_inp_data << "Number of variables: " << *ndof << "\n";
    return;}
  else if (m_model == "Cneq" || m_model == "TCneq") {
 
   file.open(getInputFiles().c_str());
std::cout << "aperto file" << std::endl;
   setMixtureFileName();
  std::cout << "set mixture file name" << std::endl;
   setnbSpecies();
std::cout << "num specie" << std::endl;
 
   setSize();
std::cout << "set size" << std::endl;
   setGlobalIndex();
std::cout << "set global index" << std::endl;
   setnbDof();
std::cout << "set ndof" << std::endl;

   setSpeciesNames();
std::cout << "set nomi di specie" << std::endl;
   setMolWeights();
std::cout << "set peso molecolare" << std::endl;
   setFormEnthalp();
   setVibrTemp();
   setSpecHeatRatio();
   setMolTypes();

   file.close();}
  else {
  //re_inp_data << "Model not implemented\n ";
   exit(1);
  }
}

//--------------------------------------------------------------------------//

std::string ChemInfo::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  std::cout << m_inputFile[0] << std::endl;
  return name;
}

//--------------------------------------------------------------------------//

void ChemInfo::setModel()
{
  model->at(0) = m_model;
 // re_inp_data << "Thermodynamic model: " << m_model << "\n";
}

//--------------------------------------------------------------------------//

void ChemInfo::setMixture() 
{
  mixture->at(0) = m_mixture;
 // re_inp_data << "File mixture data: " << getInputFile() << "\n";
}

//--------------------------------------------------------------------------//

void ChemInfo::setMixtureFileName()
{
  std::cout << "dentro a mixturefile name" << std::endl;
  file >> mixtureFile_name;
  std::cout << mixtureFile_name << std::endl;
  if (mixtureFile_name != m_mixture){
  // re_inp_data << "error data file is wrong\n";
 //  re_inp_data << mixtureFile_name << " " << m_mixture << "\n";
   exit(1);}
}

//--------------------------------------------------------------------------//

void ChemInfo::setnbSpecies()
{
  file >> *nsp;
  std::cout << *nsp << std::endl; 
  // re_inp_data << "Number of the species: " << *nsp << "\n";
}

//--------------------------------------------------------------------------//

void ChemInfo::setSpeciesNames()
{
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> name->at(ISP);
  // re_inp_data << "Chemical species: " << name->at(ISP) << ", ";
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setMolWeights()
{
 std::cout << "dentro a mol weight" << std::endl;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> mm->at(ISP);
  // re_inp_data << "MM: " << mm->at(ISP) << ", ";
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setFormEnthalp() 
{
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> hf->at(ISP);
  // re_inp_data << "HF: " << hf->at(ISP) << ", ";
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setVibrTemp()
{
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> thev->at(ISP);
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setSpecHeatRatio()
{
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> gams->at(ISP);
  // re_inp_data << "GAMS: " << gams->at(ISP) << ", ";
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setMolTypes()
{
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> typemol->at(ISP);
  }
}

//--------------------------------------------------------------------------//

void ChemInfo::setSize()
{
  std::cout << *nsp << std::endl;
  model->resize(*nsp);
  mixture->resize(*nsp);
  name->resize(*nsp);
  mm->resize(*nsp);
  hf->resize(*nsp);
  thev->resize(*nsp);
  gams->resize(*nsp);
  typemol->resize(*nsp);
} 

//--------------------------------------------------------------------------//

void ChemInfo::setGlobalIndex()
{
  *ie = m_IE + (*nsp); // = NSP+1
  *ix = m_IX + (*nsp); // = NSP+2
  *iy = m_IY + (*nsp); // = NSP+3
  *iev = m_IEV + (*nsp); // = NSP+4
}

//--------------------------------------------------------------------------//

void ChemInfo::setnbDof()
{
  if (m_model=="Cneq") { (*ndof)=(*nsp)+(*ndim)+1; }
  else if (m_model =="TCneq") { (*ndof)=(*nsp)+(*ndim)+2; }
  else {
  // re_inp_data << "Model not implemented\n";
   exit(1);
  }
  //re_inp_data << "Number of variables: " << *ndof << "\n";
}

//--------------------------------------------------------------------------//

void ChemInfo::setPhysicsData()
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
  name = PhysicsData::getInstance().getData <std::vector<std::string> > ("NAME");
  mm = PhysicsData::getInstance().getData <std::vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <std::vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <std::vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <std::vector<double> > ("GAMS");
  typemol = PhysicsData::getInstance().getData <std::vector<std::string> > ("TYPEMOL");
}
//--------------------------------------------------------------------------//

} //namespace ShockFitting
