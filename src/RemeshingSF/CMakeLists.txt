LIST ( APPEND RemeshingSF_files
BndryNodePtr.cxx
BndryNodePtr.hh
ChangeBndryPtr.cxx
ChangeBndryPtr.hh
CoNorm.cxx
CoNorm.hh
CoNorm4B.cxx
CoNorm4B.hh
CoNormPG.cxx
CoNormPG.hh
CoNorm4Ar.cxx
CoNorm4Ar.hh
CoNorm4TCneq.cxx
CoNorm4TCneq.cxx
DummyRemeshing.cxx
DummyRemeshing.hh
FndPhPs.cxx
FndPhPs.hh
RdDpsEq.cxx
RdDpsEq.hh
Rdshp.cxx
Rdshp.hh
ShpDpndnc.cxx
ShpDpndnc.hh
)

LIST ( APPEND RemeshingSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( RemeshingSF )

#SF_WARN_ORPHAN_FILES()

