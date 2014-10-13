LIST ( APPEND RemeshingSF_files
BndryNodePtr.hh
BndryNodePtr.cxx
#FndPhPs.hh
#FndPhPs.cxx
RdDpsEq.hh
RdDpsEq.cxx
)

LIST ( APPEND RemeshingSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( RemeshingSF )

#SF_WARN_ORPHAN_FILES()

