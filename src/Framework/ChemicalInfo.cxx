// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ChemicalInfo.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//


ChemicalInfo::ChemicalInfo(const std::string& objectName) : 
  Counter(),
  ConfigObject(objectName)
{
  m_Qref = 1e0;
  addOption("Qref",&m_Qref,
            "reference speed");
  m_IE = 0;
  addOption("IE",&m_IE,
             "Global index");
  m_IX = 1;
  addOption("IX",&m_IX,
             "Global index");
  m_IY = 2;
  addOption("IY",&m_IY,
             "Global index");
  m_IEV = 3;
  addOption("IEV",&m_IEV,
             "Global index");
  m_model = "dummyModel";
  addOption("Model",&m_model,
	    "Gas model");
  m_mixture = 1;
  addOption("Mixture",&m_mixture,
	    "Gas Mixture");
  m_inputFile = vector<string>();
  addOption("InputFile",&m_inputFile,
            "Gas mixture data input file");
}

//--------------------------------------------------------------------------//

string ChemicalInfo::m_model="dummymodel";

//--------------------------------------------------------------------------//

string ChemicalInfo::m_mixture="dummymixture";

//--------------------------------------------------------------------------//

double ChemicalInfo::m_Qref=0;

//--------------------------------------------------------------------------//

ChemicalInfo::~ChemicalInfo()
{
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setup()
{
  LogToScreen(VERBOSE, "ChemicalInfo::setup() => start\n");

  LogToScreen(VERBOSE, "ChemicalInfo::setup() => end\n");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::unsetup()
{
  LogToScreen(VERBOSE, "ChemicalInfo::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::read()
{
  LogToScreen(INFO, "ChemicalInfo::read()\n");

  setPhysicsData();

  logfile.Open(getClassName());

  logfile("Thermodynamic model: ",m_model, "\n");

  if (m_model == "PG") {
   (*nsp) = 1; (*ndof)=4;
   logfile("Number of variables: ",*ndof,"\n");
    return;}

  else if (m_model == "Cneq" || m_model == "TCneq") {

   logfile("File mixture data: ",getInputFiles(),"\n");

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
   setnbMol();

   file.close();
  }

  else {
   logfile("Model not implemented");
   exit(1);
  }

  logfile.Close();

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

void ChemicalInfo::setMixtureFileName()
{
  file >> dummystr;
  file >> mixtureFile_name;
  if (mixtureFile_name != m_mixture){
   logfile("ChemicalInfo::error => data file is wrong");
   logfile(mixtureFile_name, m_mixture.c_str(),"\n");
   exit(1);
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setnbSpecies()
{
  file >> dummystr;
  file >> *nsp;
  logfile( "Number of the species: ",*nsp,"\n");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setSpeciesNames()
{
  file >> dummystr;
  logfile("Chemical species: ");
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> name->at(ISP);
   logfile(name->at(ISP)," ");
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setMolWeights()
{
  file >> dummystr;
  logfile("\nMM: ");
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> mm->at(ISP);
   logfile(mm->at(ISP)," ");
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setFormEnthalp() 
{
  file >> dummystr;
  logfile("\nHF: ");
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> hf->at(ISP);
   logfile(hf->at(ISP)," ");
   if (m_mixture=="ar4") { hf->at(ISP) = hf->at(ISP)/(m_Qref*m_Qref);}
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
  logfile("\nGAMS:");
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   file >> gams->at(ISP);
   logfile(gams->at(ISP)," ");
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

void ChemicalInfo::setnbMol()
{
  *nmol = 0;
  for(unsigned ISP=0; ISP<(*nsp); ISP++) {
   if(typemol->at(ISP)=="B") { ++(*nmol); }
  }
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setSize()
{
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
  if (m_model=="Cneq") 
   { (*ndof)=(*nsp)+PhysicsInfo::getnbDim()+1; }
  else if (m_model =="TCneq") 
   { (*ndof)=(*nsp)+PhysicsInfo::getnbDim()+2; }
  else {
   logfile("Model not implemented");
   exit(1);
  }
  logfile("Number of variables: ",*ndof,"\n");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  ie = PhysicsData::getInstance().getData <unsigned> ("IE");
  ix = PhysicsData::getInstance().getData <unsigned> ("IX");
  iy = PhysicsData::getInstance().getData <unsigned> ("IY");
  iev = PhysicsData::getInstance().getData <unsigned> ("IEV");
  nmol = PhysicsData::getInstance().getData <unsigned> ("NMOL");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  name = 
    PhysicsData::getInstance().getData <vector<string> > ("NAMESP");
  mm = PhysicsData::getInstance().getData <vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <vector<double> > ("THEV");  
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  typemol =
    PhysicsData::getInstance().getData <vector<string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

void ChemicalInfo::configure(SConfig::OptionMap& cmap, 
			     const std::string& prefix)
{
  LogToScreen(VERBOSE, "ChemicalInfo::configure() => start\n");
  SConfig::ConfigObject::configure(cmap, prefix);
  LogToScreen(VERBOSE, "ChemicalInfo::configure() => end\n");
}  
  
//--------------------------------------------------------------------------//

} //namespace ShockFitting
