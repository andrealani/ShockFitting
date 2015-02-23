// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/ShockFileConverter.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"
#include "SConfig/Factory.hh"
#include "SConfig/ObjectProvider.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ShockFileConverter, Converter> 
shockFileConverterProv("ShockFileConverter");

//----------------------------------------------------------------------------//

ShockFileConverter::ShockFileConverter(const std::string& objectName) :
 Converter(objectName)
{
  m_shInputfile = "DummyshFile";
  addOption("InputFile",&m_shInputfile,
            "File containing the shock polyline and the downstream state");
  m_nbDof = 0;
  addOption("nbDof",&m_nbDof,
            "Number of variables inside the shock input file");
  m_nbShocks = 0;
  addOption("nbShocks",&m_nbShocks,
            "Number of shocks");
  m_nbSpecPoints = 0;
  addOption("nbSpecPoints",&m_nbSpecPoints,
            "Number of special points");
  m_typeSpecialPoint = "DymmySpecialPoint";
  addOption("TypeSpecPoints",&m_typeSpecialPoint,
            "Type of the special points");
  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

ShockFileConverter::~ShockFileConverter()
{
}

//----------------------------------------------------------------------------//

void ShockFileConverter::setup()
{
  LogToScreen(VERBOSE, "ShockFileConverter::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "ShockFileConverter::setup() => end\n");
}

//----------------------------------------------------------------------------//

void ShockFileConverter::unsetup()
{
  LogToScreen(VERBOSE, "ShockFileConverter::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void ShockFileConverter::configure(OptionMap& cmap, const std::string& prefix)
{
  Converter::configure(cmap, prefix);

  // assign strings on input.case file to variable transformer object
  m_prim2param.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;

  if (ConfigFileReader::isFirstConfig()) {
   m_prim2param.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                             getProvider(m_prim2param.name())
                             ->create(m_prim2param.name()));
  }

  // configure variable transformer object
  configureDeps(cmap, m_prim2param.ptr().get());
}

//----------------------------------------------------------------------------//

void ShockFileConverter::convert()
{
  LogToScreen(INFO, "ShockFileConverter::convert()\n");

  // dummy string reading first lines of the shock input file
  string dummystring;

  // dummy variables reading digital values
  double dummyvar;

  m_prim.resize(m_nbDof);
  ZRoeShu.resize(m_nbDof);
  ZRoeShd.resize(m_nbDof);
  XYSh.resize(PhysicsInfo::getnbDim());

  // transform upstream variables
  if(m_modelTransf == "Pg") { 
   m_prim.at(0) = ReferenceInfo::getpref();
   // the free-stream speed is concordant with the x-axis
   if(ReferenceInfo::speedDirectionXaxis()) { 
    m_prim.at(1) = ReferenceInfo::geturef(); }
   // the free-stream speed is pointed as -x
   else { m_prim.at(1) = -ReferenceInfo::geturef(); }
   m_prim.at(2) = 0;
   m_prim.at(3) = ReferenceInfo::getTref();
  }
  else if (m_modelTransf == "TCneq") { 
   for(unsigned i=0; i<ReferenceInfo::getrhoiref().size(); i++) {
     m_prim.at(i) = ReferenceInfo::getrhoiref().at(i);
   }
   // the free-stream speed is concordant with the x-axis
   if(ReferenceInfo::speedDirectionXaxis()) {
    m_prim.at(ReferenceInfo::getrhoiref().size()+0) =
     ReferenceInfo::geturef();
   }
   // the free-stream speed is pointed as -x
   else { 
    m_prim.at(ReferenceInfo::getrhoiref().size()+0) =
     -ReferenceInfo::geturef();
   }
   m_prim.at(ReferenceInfo::getrhoiref().size()+1) = 0.0;
   m_prim.at(ReferenceInfo::getrhoiref().size()+2) = ReferenceInfo::getTref();
   m_prim.at(ReferenceInfo::getrhoiref().size()+3) = ReferenceInfo::getTref();
  }

  // command object transforming upstream variables
  // actually the shock points coordinates are not modified
  m_prim2param.ptr()->transform(&m_prim,&XYSh,&ZRoeShu);

  // transform downstream variables
  // compute the number of shock points
  shockdat.open(m_shInputfile.c_str());

  if(!shockdat) { 
   cout << "ShockdatCreator::error => shock input file does not exist\n";
   cout << "                          try to check its name in the input.case\n";
   exit(1);
  }

  while(dummystring!="ZONE") { shockdat >> dummystring; }

  shockdat >> dummystring; // read DT = (DOUBLE DOUBLE ..)

  nbShockPoints = 0;
  while (!shockdat.eof()) {
   for(unsigned i=0; i<PhysicsInfo::getnbDim(); i++) { shockdat >> dummyvar; }
   for(unsigned i=0; i<m_nbDof; i++) { shockdat >> dummyvar; }
   ++nbShockPoints;
  }
  --nbShockPoints; 

  shockdat.close();
  
  shockdat.clear();
  dummystring = "dummy";

  // read downstream conditions
  shockdat.open(m_shInputfile.c_str());

  while(dummystring!="ZONE") { shockdat >> dummystring; }

  shockdat >> dummystring; // read DT = (DOUBLE DOUBLE ..) 

  shfile = fopen("sh00.dat","w");
  fprintf(shfile,"%1u %s",m_nbShocks,"\n");
  fprintf(shfile,"%u %s",nbShockPoints," S\n");

  for(unsigned IPOIN=0; IPOIN<nbShockPoints; IPOIN++) {
   for(unsigned IV=0; IV<PhysicsInfo::getnbDim(); IV++) { shockdat >> XYSh.at(IV); }
   for(unsigned IV=0; IV<m_nbDof; IV++) { shockdat >> m_prim.at(IV); }

   if((IPOIN == 0) || (IPOIN == nbShockPoints-1)) {
     // set XYSh(0) to 0 for the OPY special point
     if(m_typeSpecialPoint == "OPY") {XYSh.at(0) = 0;}
   }

   // command object transforming downstream variables
   m_prim2param.ptr()->transform(&m_prim,&XYSh,&ZRoeShd);

   for(unsigned IV=0; IV<PhysicsInfo::getnbDim(); IV++) {
    fprintf(shfile,"%20.15F %s",XYSh.at(IV)," ");}
   for(unsigned IV=0; IV<m_nbDof; IV++) {
    fprintf(shfile,"%20.15F %s", ZRoeShd.at(IV)," ");}
   for(unsigned IV=0; IV<m_nbDof; IV++) {
    fprintf(shfile,"%20.15F %s", ZRoeShu.at(IV)," ");}    
   fprintf(shfile,"%s","\n");
  }

  fprintf(shfile,"%1u %s",m_nbSpecPoints,"\n");
  for(unsigned i=0; i<m_nbSpecPoints; i++) {
   if     (m_typeSpecialPoint == "OPY") 
    { fprintf(shfile,"%s","OPY\n"); } 
   else { cout << "ShockdatCreator::warning => it is not implemented for this special\n";
          cout << "                            special point yet\n";
          exit(1); }
   // this printing has to be modified for nbSpecPoints>1
   fprintf(shfile,"%s %u %s"," 1", i+1,"\n");
  }

  fclose(shfile);  
  shockdat.close();
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

