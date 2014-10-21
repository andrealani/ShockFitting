LIST ( APPEND WritingMeshSF_files
WriteTriangle.cxx
WriteTriangle.hh
)

LIST ( APPEND WritingMeshSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( WritingMeshSF )

#SF_WARN_ORPHAN_FILES()
