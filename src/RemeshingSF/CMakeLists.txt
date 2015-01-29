LIST ( APPEND RemeshingSF_files
BndryFacePtrFreez.cxx
BndryFacePtrFreez.hh
BndryNodePtr.cxx
BndryNodePtrFreez.cxx
BndryNodePtr.hh
BndryNodePtrFreez.hh
ChangeBndryPtr.cxx
ChangeBndryPtr.hh
CoIntrPnt.cxx
CoIntrPnt.hh
CoNorm.cxx
CoNorm.hh
CoNorm4Ar.cxx
CoNorm4Ar.hh
CoNorm4Pg.cxx
CoNorm4Pg.hh
CoNorm4TCneq.cxx
CoNorm4TCneq.hh
CoPntDispl.cxx
CoPntDispl.hh
DummyRemeshing.cxx
DummyRemeshing.hh
FindBEdg.cxx
FindBEdg.hh
FndPhPs.cxx
FndPhPs.hh
FixMshSps.cxx
FixMshSps.hh
RdDps.cxx
RdDps.hh
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

