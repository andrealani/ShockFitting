LIST ( APPEND TestAPISF_files testAPI.cxx )

# AL: the following line assumes that the module dirs (SF_MODULES_LIST) share the same name as the corresponding libs 
LIST ( APPEND TestAPISF_libs ${SF_KERNEL_LIBS} ${SF_KERNEL_STATIC_LIBS}
FieldInterpolatorSF
FileProcessingSF
VariableTransformerSF
)

SF_ADD_PLUGIN_APP ( TestAPISF )
 
###############################################################################

#SF_WARN_ORPHAN_FILES()

