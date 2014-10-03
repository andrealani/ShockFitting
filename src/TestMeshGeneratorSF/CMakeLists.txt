LIST ( APPEND TestMeshGeneratorSF_files testMeshGenerator.cxx )

# AL: the following line assumes that the module dirs (SF_MODULES_LIST) share the same name as the corresponding libs 
LIST ( APPEND TestMeshGeneratorSF_libs ${SF_KERNEL_LIBS} ${SF_KERNEL_STATIC_LIBS}
MeshGeneratorSF
)

SF_ADD_PLUGIN_APP ( TestMeshGeneratorSF )

###############################################################################

#SF_WARN_ORPHAN_FILES()
