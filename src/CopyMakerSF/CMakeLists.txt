LIST ( APPEND CopyMakerSF_files
DummyCopyMaker.cxx
DummyCopyMaker.hh
CopyRoeValues1.cxx
CopyRoeValues1.hh
CopyRoeValues2.cxx
CopyRoeValues2.hh
MeshBackup.cxx
MeshBackup.hh
MeshRestoring.cxx
MeshRestoring.hh
)

LIST ( APPEND CopyMakerSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( CopyMakerSF )

#SF_WARN_ORPHAN_FILES()

