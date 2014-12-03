LIST ( APPEND MeshGeneratorSF_files
DummyMeshGenerator.cxx
DummyMeshGenerator.hh
ReadTriangle.cxx
ReadTriangle.hh
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

#SF_WARN_ORPHAN_FILES()


