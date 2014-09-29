# file for general configuration for the whole framework
#
# config options that affect all the framework libraries and components
# must be handled in this directory

project(config)

ADD_SUBDIRECTORY ( sysinfo )

string ( TOUPPER "${CMAKE_BUILD_TYPE}" SF_BUILD_TYPE )

#ShockFitting_config.h
configure_file(${config_SOURCE_DIR}/ShockFitting_config.h.in ${ShockFitting_BINARY_DIR}/${SF_CONFIG_HH_FILENAME})

# create the header with the svn version
configure_file(${config_SOURCE_DIR}/ShockFitting_svnversion.hh.in ${ShockFitting_BINARY_DIR}/${SF_SVN_HH_FILENAME})
