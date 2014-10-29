LIST ( APPEND VariableTransformerSF_files
ComputeTv.cxx
ComputeTv.hh
DummyVariableTransformer.cxx
DummyVariableTransformer.hh
Param2Prim.cxx
Param2Prim.hh
Param2PrimPgAdimensional.cxx
Param2PrimPgAdimensional.hh
Param2PrimPgDimensional.cxx
Param2PrimPgDimensional.hh
Param2PrimTCneqDimensional.cxx
Param2PrimTCneqDimensional.hh
Prim2Param.cxx
Prim2Param.hh
Prim2ParamPgAdimensional.cxx
Prim2ParamPgAdimensional.hh
Prim2ParamPgDimensional.cxx
Prim2ParamPgDimensional.hh
Prim2ParamTCneqDimensional.cxx
Prim2ParamTCneqDimensional.hh
VibrEnergy.cxx
VibrEnergy.hh
)

LIST ( APPEND VariableTransformerSF_libs Framework SConfig )

SF_ADD_PLUGIN_LIBRARY ( VariableTransformerSF )

#SF_WARN_ORPHAN_FILES()

