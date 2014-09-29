### Check for Linux
IF ( UNIX AND NOT APPLE )
  IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    SET ( SF_OS_LINUX 1 )
  ELSE()
    SET ( SF_OS_UNRECOGNIZED_REASON "Unrecognized UNIX type : ShockFitting has only been tested for Linux or MacOSX type UNIX'es")
  ENDIF()
ENDIF()

### Check for Apple MacOSX
IF ( APPLE )
  IF ( UNIX AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    SET ( SF_OS_MACOSX 1 )
  ELSE()
    SET ( SF_OS_UNRECOGNIZED_REASON "Unrecognized APPLE type : ShockFitting has only been tested  only with Apple MacOSX ( Darwin ) systems.")
  ENDIF()
ENDIF()

### Check for Windows
IF ( WIN32 )
  IF ( MSVC )
    SET ( SF_OS_WINDOWS 1 )
  ELSE()
    SET ( SF_OS_UNRECOGNIZED_REASON "Unrecognized WINDOWS type : ShockFitting has only been tested with Win32 and MSVC compiler.")
  ENDIF()
ENDIF()

### FINAL MESSAGE
IF ( SF_OS_UNRECOGNIZED_REASON )
  SET ( SF_OS_UNRECOGNIZED 1 )
  IF ( NOT SF_OS_SKIP_TEST )
    SET ( FULL_MSG "${SF_OS_UNRECOGNIZED_REASON} Set CMake variable SF_SKIP_OS_TEST to avoid this error" )
    MESSAGE ( FATAL_ERROR ${FULL_MSG} )
  ELSE()
    MESSAGE ( STATUS ${SF_OS_UNRECOGNIZED_REASON} )
    MESSAGE ( STATUS "Nevertheless we try to continue ..." )
  ENDIF()
ENDIF()
