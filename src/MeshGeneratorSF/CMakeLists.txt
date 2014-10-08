LIST ( APPEND MeshGeneratorSF_files
ChemicalInfo.hh
ChemicalInfo.cxx
MeshInfo.hh
MeshInfo.cxx
PhysicsInfo.hh
PhysicsInfo.cxx
ReadTriangle.hh
ReadTriangle.cxx
ReferenceInfo.hh
ReferenceInfo.cxx
ReSdwInfo.hh
ReSdwInfo.cxx
)

LIST ( APPEND MeshGeneratorSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( MeshGeneratorSF )

#SF_WARN_ORPHAN_FILES()


