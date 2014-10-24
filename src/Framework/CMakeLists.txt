LIST ( APPEND Framework_files
BaseShockFitting.cxx
BaseShockFitting.hh
ChemicalConsts.hh
ChemicalInfo.cxx
ChemicalInfo.hh
Clist.hh
Compatibility.hh
#ComputeStateDpsAdmin.cxx
#ComputeStateDpsAdmin.hh
Connectivity.hh
FieldInterpolator.cxx
FieldInterpolator.hh
FileLogManip.hh
FileProcessing.cxx
FileProcessing.hh
IOFunctions.hh
Log.hh
MeshData.hh
MeshGenerator.cxx
MeshGenerator.hh
PhysicsData.hh
PhysicsInfo.cxx
PhysicsInfo.hh
ReadInterpolateWrite.cxx
ReadInterpolateWrite.hh
ReferenceInfo.cxx
ReferenceInfo.hh
Remeshing.cxx
Remeshing.hh
ShockFittingAPI.cxx
ShockFittingAPI.hh
ShockFittingObj.cxx
ShockFittingObj.hh
ShockFittingManager.cxx
ShockFittingManager.hh
StandardShockFitting.hh
StandardShockFitting.cxx
#StateUpdater.cxx
#StateUpdater.hh
VariableTransformer.cxx
VariableTransformer.hh
WritingMesh.cxx
WritingMesh.hh
)

LIST ( APPEND Framework_libs SConfig MathTools)

SF_ADD_KERNEL_LIBRARY ( Framework )

#SF_WARN_ORPHAN_FILES()

