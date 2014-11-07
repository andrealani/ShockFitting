LIST ( APPEND TestStandardSF_files testStandardSF.cxx )

# AL: the following line assumes that the module dirs (SF_MODULES_LIST) share the same name as the corresponding libs 
LIST ( APPEND TestStandardSF_libs ${SF_KERNEL_LIBS} ${SF_KERNEL_STATIC_LIBS}
MeshGeneratorSF
RemeshingSF
WritingMeshSF
CFDSolverSF
ConverterSF
CopyMakerSF
StateUpdaterSF
VariableTransformerSF
)

SF_ADD_PLUGIN_APP ( TestStandardSF )

###############################################################################

#SF_WARN_ORPHAN_FILES()
