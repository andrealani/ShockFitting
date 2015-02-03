LIST ( APPEND ConverterSF_files
CFmesh2StartingTriangle.cxx
CFmesh2StartingTriangleFreez.cxx
CFmesh2StartingTriangle.hh
CFmesh2StartingTriangleFreez.hh
CFmesh2Triangle.cxx
CFmesh2Triangle.hh
CFmesh2TriangleFreez.cxx
CFmesh2TriangleFreez.hh
DummyConverter.cxx
DummyConverter.hh
ShockFileConverter.cxx
ShockFileConverter.hh
Triangle2CFmesh.cxx
Triangle2CFmeshFreez.cxx
Triangle2CFmeshFreez.hh
Triangle2CFmesh.hh
)

LIST ( APPEND ConverterSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( ConverterSF )

#SF_WARN_ORPHAN_FILES()

