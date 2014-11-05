LIST ( APPEND StateUpdaterSF_files
DummyStateUpdater.cxx
DummyStateUpdater.hh
CoDc.cxx
CoDc.hh
ComputeStateDps.cxx
ComputeStateDps.hh
ComputeStateDpsAr.cxx
ComputeStateDpsAr.hh
ComputeStateDpsPg.cxx
ComputeStateDpsPg.hh
ComputeStateDpsTCneq.cxx
ComputeStateDpsTCneq.hh
CopyRoeValues.cxx
CopyRoeValues.hh
CopyRoeValues1_0.cxx
CopyRoeValues1_0.hh
CoShock.cxx
CoShock.hh
FixStateSps.cxx
FixStateSps.hh
Interp.cxx
Interp.hh
MoveDps.cxx
MoveDps.hh
MoveDpsPg.cxx
MoveDpsPg.hh
MoveDpsAr.cxx
MoveDpsAr.hh
MoveDpsTCneq.cxx
MoveDpsTCneq.hh
)

LIST ( APPEND StateUpdaterSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( StateUpdaterSF )

#SF_WARN_ORPHAN_FILES()
