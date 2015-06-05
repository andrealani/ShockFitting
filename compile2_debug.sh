#!/bin/bash
cmake .. -DMPI_HOME=/software/alternate/coolfluid/cf2/2014.11/mpich2 -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_INSTALL_PREFIX=$PWD/../INSTALL -DSF_BUILD_Framework_API=ON
