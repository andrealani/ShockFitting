LIST ( APPEND ConverterSF_files
CFmesh2StartingTriangle.cxx
CFmesh2StartingTriangle.hh
CFmesh2Triangle.cxx
CFmesh2Triangle.hh
DummyConverter.cxx
DummyConverter.hh
Triangle2CFmesh.cxx
Triangle2CFmesh.hh
ShockFileConverter.cxx
ShockFileConverter.hh
)

LIST ( APPEND ConverterSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( ConverterSF )

#SF_WARN_ORPHAN_FILES()

