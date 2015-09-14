LIST ( APPEND ShockDetectorSF_files
ComputeGradient.cxx
ComputeGradient.hh
ComputeShockPointsState.cxx
ComputeShockPointsState.hh
DensityGradientMaximaPg2D.cxx
DensityGradientMaximaPg2D.hh
DetectorAlgorithm.cxx
DetectorAlgorithm.hh
DummyShockDetector.cxx
DummyShockDetector.hh
FindShockSurfaceNormals2D.cxx
FindShockSurfaceNormals2D.hh
GnoffoShockSensor.cxx
GnoffoShockSensor.hh
IdentifySpecialPoints.cxx
IdentifySpecialPoints.hh
NormalMachNumberPg2D.cxx
NormalMachNumberPg2D.hh
RedistributeShockPoints.cxx
RedistributeShockPoints.hh
)

LIST ( APPEND ShockDetectorSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( ShockDetectorSF )

#SF_WARN_ORPHAN_FILES()
