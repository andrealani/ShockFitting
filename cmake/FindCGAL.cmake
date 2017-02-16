# - Find CGAL
# This module looks for CGAL (Computational Geometry Algorithms Library) support
# Sets:
# CGAL_INCLUDE_DIR  = where CGAL headers can be found
# CGAL_LIBRARIES    = the list of library to link against
# SF_HAVE_CGAL      = set to true after finding the library
#

SET_TRIAL_INCLUDE_PATH ("") # clear include search path
SET_TRIAL_LIBRARY_PATH ("") # clear library search path

# try in user defined paths first
ADD_TRIAL_INCLUDE_PATH( ${CGAL_HOME}/include )
ADD_TRIAL_INCLUDE_PATH( $ENV{CGAL_HOME}/include )
ADD_TRIAL_INCLUDE_PATH( $ENV{CGAL_HOME}/include/CGAL )

FIND_PATH(CGAL_INCLUDE_DIR
          NAMES cartesian_homogeneous_conversion.h
          PATH_SUFFIXES CGAL
          PATHS
          ${TRIAL_INCLUDE_PATHS}  NO_DEFAULT_PATH )

# try in these paths first and then the system ones
IF ( NOT CGAL_INCLUDE_DIR )
  FIND_PATH(CGAL_INCLUDE_DIR
            NAMES cartesian_homogeneous_conversion.h
            PATH_SUFFIXES CGAL
            PATHS
            /usr/local/include
            /usr/include
	    /opt/local/include )
ENDIF()

ADD_TRIAL_LIBRARY_PATH(${CGAL_HOME}/lib )
ADD_TRIAL_LIBRARY_PATH($ENV{CGAL_HOME}/lib )

FIND_LIBRARY(CGAL_LIBRARY
             NAMES CGAL
             PATH_SUFFIXES cgal/lib CGAL/lib
             PATHS ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH )

FIND_LIBRARY(CGAL_Core_LIBRARY
             NAMES CGAL_Core
             PATH_SUFFIXES cgal/lib CGAL/lib
             PATHS ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH )

FIND_LIBRARY(CGAL_ImageIO_LIBRARY
             NAMES CGAL_ImageIO
             PATH_SUFFIXES cgal/lib CGAL/lib
             PATHS ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH )

ADD_DEFINITIONS  ( -DSF_HAVE_CGAL )
IF ( CGAL_INCLUDE_DIR AND CGAL_LIBRARY )
  SET(SF_HAVE_CGAL 1 CACHE BOOL "Found CGAL library")
ELSE()
  SET(SF_HAVE_CGAL 0 CACHE BOOL "Not found CGAL library")
ENDIF()

IF ( CGAL_LIBRARY )
  LIST ( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} )
ENDIF()
IF ( CGAL_Core_LIBRARY )
  LIST ( APPEND CGAL_LIBRARIES ${CGAL_Core_LIBRARY} )
ENDIF()
IF ( CGAL_ImageIO_LIBRARY )
  LIST ( APPEND CGAL_LIBRARIES ${CGAL_ImageIO_LIBRARY} )
ENDIF() 

MARK_AS_ADVANCED(
  CGAL_INCLUDE_DIR
  CGAL_LIBRARIES
  SF_HAVE_CGAL
)

IF ( SF_HAVE_CGAL )
  INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
ENDIF()

LOG ( "  SF_HAVE_CGAL        : [${SF_HAVE_CGAL}]" )
LOG ( "  CGAL_INCLUDE_DIR    : [${CGAL_INCLUDE_DIR}]" )
LOG ( "  CGAL_LIBRARIES      : [${CGAL_LIBRARIES}]" )
