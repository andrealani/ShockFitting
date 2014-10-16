LIST ( APPEND Framework_files
BaseShockFitting.cxx
BaseShockFitting.hh
ChemicalConsts.hh
ChemicalInfo.cxx
ChemicalInfo.hh
Clist.hh
Compatibility.hh
Connectivity.hh
FieldInterpolator.cxx
FieldInterpolator.hh
FileLogManip.hh
FileProcessing.cxx
FileProcessing.hh
IOFunctions.hh
Log.hh
MeshData.hh
MeshGenerator.hh
MeshGenerator.cxx
NormalUnitVect.cxx
NormalUnitVect.hh
PhysicsData.hh
PhysicsInfo.cxx
PhysicsInfo.hh
ReadInterpolateWrite.cxx
ReadInterpolateWrite.hh
ReferenceInfo.cxx
ReferenceInfo.hh
Remeshing.hh
Remeshing.cxx
ShockFittingAPI.cxx
ShockFittingAPI.hh
ShockFittingObj.cxx
ShockFittingObj.hh
ShockFittingManager.cxx
ShockFittingManager.hh
StandardShockFitting.hh
StandardShockFitting.cxx
VariableTransformer.cxx
VariableTransformer.hh
)

LIST ( APPEND Framework_libs SConfig MathTools)

SF_ADD_KERNEL_LIBRARY ( Framework )

#SF_WARN_ORPHAN_FILES()

