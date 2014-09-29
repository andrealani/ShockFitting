LIST ( APPEND Framework_files
BaseShockFitting.cxx
BaseShockFitting.hh
Clist.hh
Compatibility.hh
Connectivity.hh
ShockFittingAPI.cxx
ShockFittingAPI.hh
ShockFittingObj.cxx
ShockFittingObj.hh
ShockFittingManager.cxx
ShockFittingManager.hh
FieldInterpolator.cxx
FieldInterpolator.hh
FileProcessing.cxx
FileProcessing.hh
IOFunctions.hh
Log.hh
ReadInterpolateWrite.cxx
ReadInterpolateWrite.hh
VariableTransformer.cxx
VariableTransformer.hh
)

LIST ( APPEND Framework_libs SConfig )

SF_ADD_KERNEL_LIBRARY ( Framework )

#SF_WARN_ORPHAN_FILES()

