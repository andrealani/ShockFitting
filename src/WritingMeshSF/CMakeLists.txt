LIST ( APPEND WritingMeshSF_files
DummyWritingMesh.cxx
DummyWritingMesh.hh
WriteBackTriangle.cxx
WriteBackTriangle.hh
WriteSdwInfo.cxx
WriteSdwInfo.hh
WriteTriangle.cxx
WriteTriangle.hh
)

LIST ( APPEND WritingMeshSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( WritingMeshSF )

#SF_WARN_ORPHAN_FILES()
