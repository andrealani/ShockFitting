# sunversion version check
FIND_PACKAGE(Subversion)

IF (Subversion_FOUND)

  Subversion_WC_INFO(${ShockFitting_SOURCE_DIR} ShockFitting)
#   MESSAGE("Current revision is ${ShockFitting_WC_REVISION}")
#   MESSAGE("svn info : ${ShockFitting_WC_INFO}")

  FIND_PROGRAM(Subversion_SVNVERSION_EXECUTABLE svnversion DOC "subversion svnversion command line client")
  MARK_AS_ADVANCED(Subversion_SVNVERSION_EXECUTABLE)

  IF(Subversion_SVNVERSION_EXECUTABLE)
  EXECUTE_PROCESS(COMMAND ${Subversion_SVNVERSION_EXECUTABLE} -n ${ShockFitting_SOURCE_DIR}
      WORKING_DIRECTORY ${ShockFitting_SOURCE_DIR}
      OUTPUT_VARIABLE ShockFitting_svnversion
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ELSE()
    LOG("Subversion svn command was found, but not svnversion.")
    SET(ShockFitting_svnversion "NOVERSION-FOUND")
  ENDIF()

ELSE (Subversion_FOUND)
  SET(ShockFitting_svnversion "NOVERSION-FOUND")
ENDIF (Subversion_FOUND)
