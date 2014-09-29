LIST ( APPEND TestIOInterpolateSF_files testIOInterpolate.cxx )

# AL: the following line assumes that the module dirs (SF_MODULES_LIST) share the same name as the corresponding libs 
LIST ( APPEND TestIOInterpolateSF_libs ${SF_KERNEL_LIBS} ${SF_KERNEL_STATIC_LIBS}
FieldInterpolatorSF
FileProcessingSF
VariableTransformerSF
)

SF_ADD_PLUGIN_APP ( TestIOInterpolateSF )
 
###############################################################################

#SF_WARN_ORPHAN_FILES()

