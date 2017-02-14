##############################################################################
# ShockFitting cmake build system
##############################################################################
# TODO
#  * packaging with the build system
#  * check precompiled headers
##############################################################################
# Command line options for cmake:
# -DCMAKE_BUILD_TYPE=DEBUG     will set the build type
# -DSF_SKIP_FORTRAN=:BOOL=ON   will not test for fortran language
# -DCMAKE_SKIP_RPATH:BOOL=ON   will remove rpath
# -G"Visual Studio 9 2008"     will generate for WIN32
# -DCMAKE_INSTALL_PREFIX=DIR   will specify where to install
##############################################################################

##############################################################################
# cmake behavior
##############################################################################

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 FATAL_ERROR )

# Tell the CMake makefile generator to not have rules depend on themselves.
# This causes extra rebuilds when the include path changes from turning a library on or off
SET(CMAKE_SKIP_RULE_DEPENDENCY 1)
# allow relaxed loop and if constructions
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

##############################################################################
# project definition and versioning
##############################################################################

PROJECT (ShockFitting C CXX )

SET (ShockFitting_VERSION_MAJOR    2014)
SET (ShockFitting_VERSION_MINOR    1)
SET (ShockFitting_VERSION "${ShockFitting_VERSION_MAJOR}.${ShockFitting_VERSION_MINOR}" )

##############################################################################
# cmake setup
##############################################################################

# disallow in-source build
IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "ShockFitting requires an out of source build.\nPlease create a separate build directory and run 'cmake path_to_ShockFitting [options]' there.")
ENDIF()

# include other configurations
SET(SF_CMAKE_DIR "${ShockFitting_SOURCE_DIR}/cmake")
SET(CMAKE_MODULE_PATH "${SF_CMAKE_DIR}" ${CMAKE_MODULE_PATH})

INCLUDE( DefineMacros )   # define macros
#INCLUDE( CTest )          # include testing support

# recreate the logfile
SET  ( PROJECT_LOG_FILE ${PROJECT_BINARY_DIR}/CMakeLogInfo.txt )
FILE ( WRITE ${PROJECT_LOG_FILE} "ShockFitting Cmake Log file\n")

MESSAGE ( STATUS "CMAKE_INSTALL_PREFIX [${CMAKE_INSTALL_PREFIX}]"  ) 

IF ( NOT DEFINED SF_INSTALL_BIN_DIR )
  SET(SF_INSTALL_BIN_DIR bin)
ENDIF ()

IF ( NOT DEFINED SF_INSTALL_LIB_DIR )
  SET(SF_INSTALL_LIB_DIR lib)
ENDIF ()

IF ( NOT DEFINED SF_INSTALL_INCLUDE_DIR )
  SET(SF_INSTALL_INCLUDE_DIR include/ShockFitting)
ENDIF ()

IF ( NOT DEFINED SF_INSTALL_DATA_DIR )
  SET(SF_INSTALL_DATA_DIR data)
ENDIF ()

# use RPATHs, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# create the dso directory for shared libraries
SET ( ShockFitting_DSO_DIR ${ShockFitting_BINARY_DIR}/dso )
FILE ( MAKE_DIRECTORY ${ShockFitting_DSO_DIR} )

# name the generated config files
SET ( SF_CONFIG_HH_FILENAME ShockFitting_config.h CACHE INTERNAL "name of configuration file" )
SET ( SF_SVN_HH_FILENAME    ShockFitting_svnversion.hh CACHE INTERNAL "name of svn version file" )

##############################################################################
# configuration options
##############################################################################

# set the default build to be OPTIM
IF ( NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE OPTIM CACHE STRING "Choose the type of build, options are: DEBUG OPTIM RELEASE." FORCE)
ENDIF()

INCLUDE ( CheckOperatingSystem  )   # check for operating system
INCLUDE ( CheckCompilerFeatures )   # check compiler features

INCLUDE ( DefineGlobalOptions   )   # add user global options
INCLUDE ( DefineBuildModes      )   # create extra build modes
INCLUDE ( DefineBuildRules      )   # defined default compilation flags and linking rules per architecture

##############################################################################
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "General Configuration ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
##############################################################################

ADD_SUBDIRECTORY ( cmake )
#ADD_SUBDIRECTORY ( doc )

##############################################################################

# our find macros

# MPI 
FIND_PACKAGE(MPI)

#include the generated config file
ADD_DEFINITIONS  ( -DSF_HAVE_CONFIG_H )
INCLUDE_DIRECTORIES ( ${ShockFitting_BINARY_DIR} )

SET(SF_HAVE_TRILIBRARY 1 CACHE BOOL "Definition of symbol for compiling Triangle as a library")
ADD_DEFINITIONS  ( -DNO_TIMER )
ADD_DEFINITIONS  ( -DANSI_DECLARATORS )
ADD_DEFINITIONS  ( -DTRILIBRARY )

# CGAL
FIND_PACKAGE(CGAL)
#IF (NOT CGAL_HOME EQUAL "")
#SET(SF_HAVE_CGAL 1 CACHE BOOL "Definition of symbol for using CGAL library")
#ENDIF()

##############################################################################
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "Kernel Configuration ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
##############################################################################

# continue to kernel modules
INCLUDE_DIRECTORIES( ${ShockFitting_SOURCE_DIR}/src )
# include in compilation path the src dir
ADD_SUBDIRECTORY ( src )

#config file generation (at the end, so everything is detected)
ADD_SUBDIRECTORY ( config )

# install the headers
INSTALL(FILES ${ShockFitting_BINARY_DIR}/${SF_CONFIG_HH_FILENAME}
              ${ShockFitting_BINARY_DIR}/${SF_SVN_HH_FILENAME}
        DESTINATION
              ${SF_INSTALL_INCLUDE_DIR})

##############################################################################
# print summary
##############################################################################

LOG ( "")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( "+ ShockFitting configuration summary ")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( " ShockFitting Release     : [${ShockFitting_VERSION}]")
LOG ( " ShockFitting SVN version : [${ShockFitting_svnversion}]")
LOG ( "---------------------------------------------------------" )
LOG ( " CMake Version         : [${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}]")
LOG ( " CMake Generator       : [${CMAKE_GENERATOR}]")
LOG ( " Build Type            : [${CMAKE_BUILD_TYPE}]")
LOG ( " Operating System      : [${CMAKE_SYSTEM}]")
LOG ( " Memory address        : [${SF_OS_BITS} bits]")
LOG ( "---------------------------------------------------------" )
LOG ( " C compiler            : [${CMAKE_C_COMPILER}]")
LOG ( " C flags               : [${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}]")
LOG ( " C++ compiler          : [${CMAKE_CXX_COMPILER}]")
LOG ( " C++ flags             : [${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}]")
LOG ( " C++ link flags        : [${CMAKE_CXX_LINK_FLAGS}]")
IF( NOT SF_SKIP_FORTRAN )
LOG ( " Fortran compiler      : [${CMAKE_Fortran_COMPILER}]")
LOG ( " Fortran flags         : [${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}]")
ENDIF()
LOG ( " Shared Linker flags   : [${CMAKE_SHARED_LINKER_FLAGS}]")

LOG ( "---------------------------------------------------------" )
LOG ( " MPI                   : [${SF_HAVE_MPI}]" )
IF(SF_HAVE_MPI)
LOG ( "    SF_MPIRUN_PROGRAM  : [${SF_MPIRUN_PROGRAM}]" )
LOG ( "    MPI_INCLUDE_DIR    : [${MPI_INCLUDE_DIR}]")
LOG ( "    MPI_LIBRARIES      : [${MPI_LIBRARIES}]")
ENDIF()
LOG ( " Kernel libs           : [${SF_KERNEL_LIBS}]")
LOG ( " Static plugin libs    : [${SF_KERNEL_STATIC_LIBS}]")
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
LOG ( "")
LOG ( "    SF_HAVE_TRILIBRARY : [${SF_HAVE_TRILIBRARY}]" )
LOG ( "    SF_HAVE_CGAL       : [${SF_HAVE_CGAL}]" )
LOG ( "    CGAL_HOME          : [${CGAL_HOME}]" )
LOG ( "    CGAL_INCLUDE_DIR   : [${CGAL_INCLUDE_DIR}]" )
LOG ( "    CGAL_LIBRARIES     : [${CGAL_LIBRARIES}]" )

# if file with orphan files exists remove it
SET ( ORPHAN_FILE "${ShockFitting_BINARY_DIR}/OrphanFiles.txt" )
IF ( EXISTS ${ORPHAN_FILE} ) 
	FILE ( REMOVE ${ORPHAN_FILE} )
ENDIF()

# if orphan files where found, put the list on the file
LIST(LENGTH SF_ORPHAN_FILES SF_LENGTH_ORPHAN_FILES)
IF (SF_LENGTH_ORPHAN_FILES)
  LOG ( " !!! WARNING !!! Orphan files where found during the configuration.")
  LOG ( " !!! WARNING !!! Check full list in file ${ORPHAN_FILE} ")
  FOREACH( AFILE ${SF_ORPHAN_FILES} )
    FILE ( APPEND ${ORPHAN_FILE} "${AFILE}\n" )
  ENDFOREACH()
ENDIF()

DUMP_VARIABLES()        # dump relevant variables in the log file
#INCLUDE(PrepareCPack)   # packaging instructions

