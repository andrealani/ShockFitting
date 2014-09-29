LIST ( APPEND FieldInterpolatorSF_files
BoxFieldInterpolator.cxx
BoxFieldInterpolator.hh
DummyFieldInterpolator.cxx
DummyFieldInterpolator.hh
)

LIST ( APPEND FieldInterpolatorSF_libs Framework SConfig )

SF_ADD_PLUGIN_LIBRARY ( FieldInterpolatorSF )

#SF_WARN_ORPHAN_FILES()

