#########################################################################################
# Generic OPTIONS
#########################################################################################

# precision real numbers
IF ( NOT SF_PRECISION_SINGLE )
OPTION(SF_PRECISION_SINGLE       "Real numbers have single precision"       OFF )
ENDIF()
IF ( NOT SF_PRECISION_LONG_DOUBLE )
OPTION(SF_PRECISION_LONG_DOUBLE   "Real numbers have long double precision" OFF )
ENDIF()
IF ( NOT SF_PRECISION_SINGLE AND NOT SF_PRECISION_LONG_DOUBLE )
OPTION(SF_PRECISION_DOUBLE   "Real numbers have double precision"           ON )
ENDIF()

# user option to search other dirs for plugins
SET ( SF_EXTRA_SEARCH_DIRS "" CACHE STRING "Full paths to extra dirs to be searched for plugin modules which maybe out of source." )

OPTION ( SF_ENABLE_MPI                "Enable MPI compilation"                  ON   )
OPTION ( SF_ENABLE_DOCS               "Enable build of documentation"           ON   )
OPTION ( SF_ENABLE_WARNINGS           "Enable lots of warnings while compiling"  ON )
OPTION ( SF_CMAKE_LIST_PLUGINS             "CMake lists the plugins"                 OFF  )

# if user disables MPI we overwrite the SF_HAVE_MPI variable
IF ( SF_ENABLE_MPI )
  IF ( SF_MPI_AVAILABLE )
    SET ( SF_HAVE_MPI 1 CACHE BOOL "User enabled MPI [FOUND]" )
  ELSE()
    SET ( SF_HAVE_MPI 0 CACHE BOOL "User enabled MPI [NOT-FOUND]" )
  ENDIF()
ELSE()
  SET ( SF_HAVE_MPI 0 CACHE BOOL "User disabled MPI" )
ENDIF ()

#########################################################################################
# PROFILING OPTIONS
#########################################################################################

# user option to add system depedent profiling
OPTION ( SF_ENABLE_PROFILING    "Enable code profiling"                 OFF )

IF(SF_ENABLE_PROFILING)

  ###########################
  # GNU gprof
  IF(SF_PROFILER_TOOL MATCHES gprof)
    IF(UNIX AND CMAKE_COMPILER_IS_GNUCC)
      SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -pg" )
      SET(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -pg" )
    ELSE(UNIX AND CMAKE_COMPILER_IS_GNUCC)
      LOG("User selected profiler [gprof] must be used with GCC compiler")
      SET( SF_PROFILER_TOOL     NOTFOUND )
    ENDIF()
  ENDIF()

  ###########################
  # google-perftools
  IF(SF_PROFILER_TOOL MATCHES google-perftools)

    FIND_PACKAGE(GooglePerftools)

    IF(SF_HAVE_GOOGLE_PERFTOOLS)
      LINK_LIBRARIES(${GOOGLE_PERFTOOLS_LIBRARIES})
    ELSE(SF_HAVE_GOOGLE_PERFTOOLS)
      LOG("User selected profiler [google-pertools] could not be found")
      SET( SF_PROFILER_TOOL     NOTFOUND )
    ENDIF()

  ENDIF()

ENDIF()

#########################################################################################
# STATIC BUILD OPTIONS
#########################################################################################

# user option to static build
OPTION ( SF_ENABLE_STATIC       "Enable static building"                OFF)

IF ( SF_ENABLE_STATIC )

  SET(BUILD_SHARED_LIBS OFF)

    IF(UNIX)
      # LINUX
      IF("${CMAKE_SYSTEM}" MATCHES Linux)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -Wl,-whole-archive <LINK_LIBRARIES> -Wl,-no-whole-archive")
      ENDIF("${CMAKE_SYSTEM}" MATCHES Linux)
      # SGI IRIX
      IF("${CMAKE_SYSTEM}" MATCHES IRIX)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -Wl,-all <LINK_LIBRARIES> -Wl,-notall")
      ENDIF("${CMAKE_SYSTEM}" MATCHES IRIX)
      # On Darwin:
      #  -all_load $convenience
      IF("${CMAKE_SYSTEM}" MATCHES Darwin)
        SET(CMAKE_CXX_LINK_EXECUTABLE
        "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> -all_load <LINK_LIBRARIES>")
      ENDIF("${CMAKE_SYSTEM}" MATCHES Darwin)
      # On Solaris 2:
      #   -z allextract $convenience -z defaultextract
    ENDIF(UNIX)

ELSE()

  SET( BUILD_SHARED_LIBS ON )

ENDIF()




