// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ConverterSF/DummyConverter.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<DummyConverter, Converter>
dummyConverterProv("DummyConverter");

//--------------------------------------------------------------------------//

DummyConverter::DummyConverter(const std::string& objectName) :
  Converter(objectName)
{
}

//--------------------------------------------------------------------------//

DummyConverter::~DummyConverter()
{

}

//--------------------------------------------------------------------------//

void DummyConverter::setup()
{
}

//--------------------------------------------------------------------------//

void DummyConverter::unsetup()
{
}

//--------------------------------------------------------------------------//

void DummyConverter::convert()
{
  std::cout << "DummyConvert::convert()\n";
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
