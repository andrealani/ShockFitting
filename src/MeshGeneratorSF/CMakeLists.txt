LIST ( APPEND MeshGeneratorSF_files
MeshInfo.hh
MeshInfo.cxx
PhysicsInfo.hh
PhysicsInfo.cxx
ReadTriangle.hh
ReadTriangle.cxx
)

LIST ( APPEND MeshGeneratorSF_libs Framework SConfig )

SF_ADD_PLUGIN_LIBRARY ( MeshGeneratorSF )

#SF_WARN_ORPHAN_FILES()


