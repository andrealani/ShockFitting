LIST ( APPEND MeshGeneratorSF_files
DummyMeshGenerator.cxx
DummyMeshGenerator.hh
MedianDualCell.cxx
MedianDualCell.hh
ReadTriangle.cxx
ReadTriangleFreez.cxx
ReadTriangle.hh
ReadTriangleFreez.hh
ReSdwInfo.cxx
ReSdwInfo.hh
)

LIST ( APPEND MeshGeneratorSF_libs Framework SConfig MathTools)

SF_ADD_PLUGIN_LIBRARY ( MeshGeneratorSF )

IF (SF_HAVE_TRILIBRARY)

LIST ( APPEND triangleLib_files
TriLibrary/triangle.c
TriLibrary/triangle.h
)
SF_ADD_PLUGIN_LIBRARY ( triangleLib )

LIST ( APPEND TriangleSF_files
TriangleExe.cxx
TriangleExe.hh
Tricall.cxx
Tricall.hh
)
LIST ( APPEND TriangleSF_libs Framework SConfig MathTools MeshGeneratorSF triangleLib )
SF_ADD_PLUGIN_LIBRARY ( TriangleSF )
ENDIF()

# example for CGAL (not necessarily to be put here)
# IF (SF_HAVE_CGAL)
# LIST ( APPEND CGALSF_includedirs ${CGAL_INCLUDE_DIR} )
# LIST ( APPEND CGALSF_libs Framework SConfig MathTools MeshGeneratorSF ${CGAL_LIBRARIES} )
# SF_ADD_PLUGIN_LIBRARY ( CGALSF )
# ENDIF()

#SF_WARN_ORPHAN_FILES()


