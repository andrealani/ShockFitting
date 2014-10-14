LIST ( APPEND RemeshingSF_files
BndryNodePtr.cxx
BndryNodePtr.hh
ChangeBndryPtr.cxx
ChangeBndryPtr.hh
CoNorm.cxx
CoNorm.hh
CoNorm4B.cxx
CoNorm4B.hh
DummyRemeshing.cxx
DummyRemeshing.hh
FndPhPs.cxx
FndPhPs.hh
RdDpsEq.cxx
RdDpsEq.hh
RdDpsEq.cxx
)

LIST ( APPEND RemeshingSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( RemeshingSF )

#SF_WARN_ORPHAN_FILES()

