LIST ( APPEND StateUpdaterSF_files
DummyStateUpdater.cxx
DummyStateUpdater.hh
CoDc.cxx
CoDc.hh
ComputeResidual.cxx
ComputeResidual.hh
ComputeStateDps4Ar.cxx
ComputeStateDps4Ar.hh
ComputeStateDps4Pg.cxx
ComputeStateDps4Pg.hh
ComputeStateDps4TCneq.cxx
ComputeStateDps4TCneq.hh
CoShock.cxx
CoShock.hh
DiscrErrorNormL1.cxx
DiscrErrorNormL1.hh
DiscrErrorNormL1weighted.cxx
DiscrErrorNormL1weighted.hh
DiscrErrorNormL2.cxx
DiscrErrorNormL2.hh
FixStateSps.cxx
FixStateSps.hh
Interp.cxx
Interp.hh
MoveDps4Ar.cxx
MoveDps4Ar.hh
MoveDps4Pg.cxx
MoveDps4Pg.hh
MoveDps4TCneq.cxx
MoveDps4TCneq.hh
)

LIST ( APPEND StateUpdaterSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( StateUpdaterSF )

#SF_WARN_ORPHAN_FILES()
