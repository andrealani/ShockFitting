LIST ( APPEND Framework_files
BaseShockFitting.cxx
BaseShockFitting.hh
BoundaryConnectivity.hh
CFDSolver.cxx
CFDSolver.hh
ChemicalConsts.hh
ChemicalInfo.cxx
ChemicalInfo.hh
Clist.hh
Compatibility.hh
Connectivity.hh
Converter.cxx
Converter.hh
CopyMaker.cxx
CopyMaker.hh
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
PhysicsConsts.hh
PhysicsData.hh
PhysicsInfo.cxx
PhysicsInfo.hh
ReadInterpolateWrite.cxx
ReadInterpolateWrite.hh
ReferenceInfo.cxx
ReferenceInfo.hh
Remeshing.cxx
Remeshing.hh
ShockDetector.cxx
ShockDetector.hh
ShockFittingAPI.cxx
ShockFittingAPI.hh
ShockFittingObj.cxx
ShockFittingObj.hh
ShockFittingManager.cxx
ShockFittingManager.hh
StandardShockFitting.cxx
StandardShockFittingBeta.cxx
StandardShockFitting.hh
StandardShockFittingBeta.hh
StateUpdater.cxx
StateUpdater.hh
VariableTransformer.cxx
VariableTransformer.hh
WritingMesh.cxx
WritingMesh.hh
)

LIST ( APPEND Framework_libs SConfig MathTools)

SF_ADD_KERNEL_LIBRARY ( Framework )

#SF_WARN_ORPHAN_FILES()

