LIST ( APPEND ConverterSF_files
CFmesh2StartingTriangle.cxx
CFmesh2StartingTriangleFreez.cxx
CFmesh2StartingTriangleFVM.cxx
CFmesh2StartingTriangle.hh
CFmesh2StartingTriangleFreez.hh
CFmesh2StartingTriangleFVM.hh
CFmesh2Triangle.cxx
CFmesh2Triangle.hh
CFmesh2TriangleFreez.cxx
CFmesh2TriangleFreez.hh
CFmesh2TriangleFVM.cxx
CFmesh2TriangleFVM.hh
DummyConverter.cxx
DummyConverter.hh
ShockFileConverter.cxx
ShockFileConverter.hh
Tecplot2StartingTriangle.cxx
Tecplot2StartingTriangle.hh
Triangle2CFmesh.cxx
Triangle2CFmeshFreez.cxx
Triangle2CFmeshFreez.hh
Triangle2CFmeshFVM.cxx
Triangle2CFmeshFVM.hh
Triangle2CFmesh.hh
Triangle2Tecplot.cxx
Triangle2Tecplot.hh
)

LIST ( APPEND ConverterSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( ConverterSF )

#SF_WARN_ORPHAN_FILES()

