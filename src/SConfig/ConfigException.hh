// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ConfigException_hh
#define ConfigException_hh  

#ifdef SF_HAVE_MPI
#include <mpi.h>
#endif

#include <exception>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This function outputs an error message only from the master process
static void printError(std::ostream& out, const std::string& error)
{
  int rank = 0;
#ifdef SF_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == 0) {
    out << error << "\n";
  }
}

// This class implements an exception indicating configuration failure
// @author Andrea Lani

class ConfigException : public std::exception {
public:
  
  // constructor
  ConfigException(const std::string& msg) :  
    std::exception(),
    m_msg(msg)
  {
  }
  
  // destructor
  virtual ~ConfigException () throw() {}
  
  // throw a char array with the error message
  const char* what() const throw() {return m_msg.c_str();}
  
private:
  
  // error message to display
  std::string m_msg;
  
};

//----------------------------------------------------------------------------//
  
} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
