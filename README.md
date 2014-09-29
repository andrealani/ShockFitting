ShockFitting
============

This will be a modular shock fitting code to be plugged to arbitrary CFD codes.


##### Installation instructions for Shock Fitting code #####

For the configuration and compilation, you will need to have cmake (version > 2.8) installed in your system.
To get started, create a directory "build" inside your CouplingTools home, then move into it:

mkdir build ; cd build

You can configure by running the command 

cmake .. -DMPI_HOME=MPIDIR -DCMAKE_CXX_COMPILER=CXX -DCMAKE_INSTALL_PREFIX=INSTALLDIR -DCF_BUILD_Framework_API=ON -DCMAKE_BUILD_TYPE=DEBUG

where:
  
CXX        : chosen C++ compiler
MPIDIR     : directory of existing MPI installation
INSTALLDIR : directory where CouplingTools libraries will be installed

The libraries can now be compiled with 

make install

Upon successful completion, you will find all shared libraries and include files (from the Framework only) respectively inside
  
INSTALLDIR/lib
INSTALLDIR/include/couplingtools/Framework

