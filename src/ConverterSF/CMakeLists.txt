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
TecplotFVM2StartingTriangle.cxx
TecplotFVM2StartingTriangle.hh
TecplotFVM2Triangle.cxx
TecplotFVM2Triangle.hh
Tecplot2StartingTriangle.cxx
Tecplot2StartingTriangle.hh
Tecplot2Triangle.cxx
Tecplot2Triangle.hh
Triangle2CFmesh.cxx
Triangle2CFmeshFreez.cxx
Triangle2CFmeshFreez.hh
Triangle2CFmesh.hh
Triangle2Tecplot.cxx
Triangle2Tecplot.hh
)

LIST ( APPEND ConverterSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( ConverterSF )

#SF_WARN_ORPHAN_FILES()

