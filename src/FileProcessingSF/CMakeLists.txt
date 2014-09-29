LIST ( APPEND FileProcessingSF_files
DummyFileProcessing.cxx
DummyFileProcessing.hh
ReadCFmesh.hh
ReadCFmesh.cxx
ReadPointCloud.hh
ReadPointCloud.cxx
WritePointCloud.hh
WritePointCloud.cxx
)

LIST ( APPEND FileProcessingSF_libs Framework SConfig )

SF_ADD_PLUGIN_LIBRARY ( FileProcessingSF )

#SF_WARN_ORPHAN_FILES()

