// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ShockFittingAPI.hh"
#include "Framework/Field.hh"
#include "SConfig/Factory.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;
using namespace ShockFitting;

//--------------------------------------------------------------------------//

void SF_create_(int* federateID)
{  
  assert(*federateID < NFEDERATES);
  obj[*federateID].reset(new ShockFittingManager());
}
  
//--------------------------------------------------------------------------//
  
void SF_configure_(int* federateID, char* input, 
		   int* argc, char*** argv)
{
  // preprocess strings (needed if coming from FORTRAN)
  string iname(input); trimRear(iname);
  
  assert(*federateID < NFEDERATES);
  obj[*federateID]->configure(iname, *argc, *argv);
  obj[*federateID]->setup();
}
  
//--------------------------------------------------------------------------//
   
void SF_process_(int* federateID)
{  
  assert(*federateID < NFEDERATES);
  obj[*federateID]->getObj().process();
}

//--------------------------------------------------------------------------//

void SF_interpolate_field_(int* federateID,
			   char* objName,
			   int* inNbElems,
			   int* inNbStates,
			   int* inStateStride,
			   int* inElementState, 
			   int* inElementStatePtr,
			   double* inState, 
			   int* outNbElems,
			   int* outNbStates,
			   int* outStateStride,
			   int* outElementState, 
			   int* outElementStatePtr,
			   double* outState)
{
  // preprocess strings (needed if coming from FORTRAN)
  string oname(objName); trimRear(oname);
  
  assert(*federateID < NFEDERATES);
  SharedPtr<ShockFittingManager> ct = obj[*federateID];
  std::vector<PAIR_TYPE(FieldInterpolator)>& list = 
    ct->getObj().getFieldInterpolatorList();
  
  for (unsigned i = 0; i < list.size(); ++i) {
    if (list[i].name() == oname) {
      // create all connectivity objects
      Connectivity inStateConn((unsigned)(*inNbElems), inElementState, inElementStatePtr);
      Connectivity outStateConn((unsigned)(*outNbElems), outElementState, outElementStatePtr);
      
      // create all field objects
      Field finState((unsigned)(*inNbStates), (unsigned)(*inStateStride), inStateConn, inState);
      Field foutState((unsigned)(*outNbStates), (unsigned)(*outStateStride), outStateConn, outState);
      
      // interpolate
      list[i].ptr()->interpolate(&finState, &foutState);
    }
  }
}

//--------------------------------------------------------------------------//
 
void SF_process_file_(int* federateID, char* objName)
{  
  // preprocess strings (needed if coming from FORTRAN)
  string oname(objName); trimRear(oname);

  assert(*federateID < NFEDERATES);
  SharedPtr<ShockFittingManager> ct = obj[*federateID];
  std::vector<PAIR_TYPE(FileProcessing)>& list = 
    ct->getObj().getFileProcessingList();
  
  for (unsigned i = 0; i < list.size(); ++i) {
    if (list[i].name() == oname) {
      list[i].ptr()->process();
    }
  }
}

//--------------------------------------------------------------------------//

void SF_destroy_(int* federateID)
{  
  obj[*federateID]->unsetup();
}

//--------------------------------------------------------------------------//
