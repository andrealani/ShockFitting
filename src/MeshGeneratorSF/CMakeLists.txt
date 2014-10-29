LIST ( APPEND MeshGeneratorSF_files
DummyMeshGenerator.cxx
DummyMeshGenerator.hh
ReadTriangle.cxx
ReadTriangle.hh
ReSdwInfo.cxx
ReSdwInfo.hh
Triangle.hh
Triangle.cxx
)

LIST ( APPEND MeshGeneratorSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( MeshGeneratorSF )

#SF_WARN_ORPHAN_FILES()


