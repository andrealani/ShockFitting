LIST ( APPEND CFDSolverSF_files
COOLFluiD.cxx
COOLFluiD.hh
)

LIST ( APPEND CFDSolverSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( CFDSolverSF )

#SF_WARN_ORPHAN_FILES()

