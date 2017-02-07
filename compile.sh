#!/bin/bash
cmake .. -DMPI_HOME=/software/alternate/coolfluid/cf2/2013.9/openmpi -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_INSTALL_PREFIX=$PWD/../INSTALL -DSF_BUILD_Framework_API=ON 
# if CGAL has to be used, add the following to the previous line with right full paths
#-DCGAL_HOME=... -DCGAL_INCLUDE_DIR=.../include -DCGAL_LIBRARIES=...dylib;...dylib;...dylib
