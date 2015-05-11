// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockFittingAPI_hh
#define ShockFitting_ShockFittingAPI_hh

//--------------------------------------------------------------------------//

#include "Framework/ShockFittingManager.hh"

//--------------------------------------------------------------------------//

/// This class defines a API for the ShockFittingObj which is usable from C++,
/// C and FORTRAN
/// 
/// @author Andrea Lani

const int NFEDERATES = 100;  
SConfig::SharedPtr<ShockFitting::ShockFittingManager> obj[NFEDERATES];
  
// --- extern definitions --- // 
#ifdef __cplusplus
extern "C" {
#endif
  
  /// create the shock fitting for the given federate
  /// @param federateID   ID of the federate (< NFEDERATES)
  void SF_create_(int* federateID);
  
  /// configure the shock fitting
  /// @param federateID   ID of the federate (< NFEDERATES)
  /// @param input       configuration filename or string
  /// @param argc        number of command line options
  /// @param argv        command line options
  void SF_configure_(int* federateID, 
		     char* input, 
		     int* argc, 
		     char*** argv);
  
  /// run the shock fitting 
  /// @param federateID   ID of the federate (< NFEDERATES)
  void SF_process_(int* federateID);
  
  /// run the shock fitting 
  /// @param federateID          ID of the federate (< NFEDERATES)
  /// @param dim                 spatial dimension (2 or 3)
  /// @param inNbBoundaries      input number of individual boundary patches
  /// @param inbInfo             input array of size @see nbBoundaries x 2, storing consecutively 
  ///                            (1) the number of elements (faces) in each boundary 
  ///                            (2) a pointer to the first element in each boundary
  /// @param inbNode             input element-node connectivity
  /// @param inbPtr              input pointer to the start of the corresponding element
  /// @param inNbElems           number of elements in the input mesh
  /// @param inNbStates          number of states in the input mesh
  /// @param inStateStride       stride for states (e.g. number of equations) in the input mesh
  /// @param inElementState      element-state connectivity, i.e. a list of all state IDs 
  ///                            for each element connectivity data for the input solution
  /// @param inElementStatePtr   array storing pointers to the beginning of each element
  ///                            for the input solution field
  /// @param inState             input solution field
  /// @param inNode              input spatial coordinates field
  /// @param outNbBoundaries     output number of individual boundary patches
  /// @param outbInfo            output array of size @see nbBoundaries x 2, storing consecutively 
  ///                            (1) the number of elements (faces) in each boundary 
  ///                            (2) a pointer to the first element in each boundary
  /// @param outbNode            output element-node connectivity
  /// @param outbPtr             output pointer to the start of the corresponding element
  /// @param outNbElems          number of elements in the output mesh
  /// @param outNbStates         number of states in the output mesh  
  /// @param outStateStride      stride for states (e.g. number of equations) in the output mesh
  /// @param outElementState     element-state connectivity, i.e. a list of all state IDs 
  ///                            for each element connectivity data for the interpolated solution
  /// @param outElementStatePtr  array storing pointers to the beginning of each element
  ///                            for the interpolated solution field
  /// @param outState            output interpolated solution field
  /// @param outNode             output spatial coordinates field
  void SF_process_field_(int* federateID, int* dim,
			 int* inNbBoundaries, int* inbInfo, int* inbNode, int* inbPtr, char* inbName,
			 int* inNbElems,
			 int* inNbStates,
			 int* inStateStride,
			 int* inElementState, 
			 int* inElementStatePtr,
			 double* inState, 
			 double* inNode, 
			 int* outNbBoundaries, int* outbInfo, int* outbNode, int* outbPtr,
			 int* outNbElems,
			 int* outNbStates,
			 int* outStateStride,
			 int* outElementState, 
			 int* outElementStatePtr,
			 double* outState,
			 double* outNode);
  
  /// Interpolate from one field and data structure into another
  /// @param federateID          ID of the federate (< NFEDERATES)
  /// @param objName             name of the actual field interpolator to use
  /// @param inNbElems           number of elements in the input mesh
  /// @param inNbStates          number of states in the input mesh
  /// @param inStateStride       stride for states (e.g. number of equations) in the input mesh
  /// @param inElementState      element-state connectivity, i.e. a list of all state IDs 
  ///                            for each element connectivity data for the input solution
  /// @param inElementStatePtr   array storing pointers to the beginning of each element
  ///                            for the input solution field
  /// @param inState             input solution field
  /// @param outNbElems          number of elements in the output mesh
  /// @param outNbStates         number of states in the output mesh  
  /// @param outStateStride      stride for states (e.g. number of equations) in the output mesh
  /// @param outElementState     element-state connectivity, i.e. a list of all state IDs 
  ///                            for each element connectivity data for the interpolated solution
  /// @param outElementStatePtr  array storing pointers to the beginning of each element
  ///                            for the interpolated solution field
  /// @param outState            interpolated solution field
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
			     double* outState);
  
  /// Process one or more input files and save the result into one or more output files
  /// @param federateID    ID of the federate (< NFEDERATES)
  /// @param objName       name of the actual file processing to use
  void SF_process_file_(int* federateID, char* objName);
  
  /// finalize the shock fitting for the given federate
  /// @param federateID    ID of the federate (< NFEDERATES)
  void SF_destroy_(int* federateID);
  
#ifdef __cplusplus
}
#endif
  
//--------------------------------------------------------------------------//

#endif // ShockFitting_ShockFittingAPI_hh
