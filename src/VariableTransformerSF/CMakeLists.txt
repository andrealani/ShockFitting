LIST ( APPEND VariableTransformerSF_files
DummyVariableTransformer.cxx
DummyVariableTransformer.hh
)

LIST ( APPEND VariableTransformerSF_libs Framework SConfig )

SF_ADD_PLUGIN_LIBRARY ( VariableTransformerSF )

#SF_WARN_ORPHAN_FILES()

