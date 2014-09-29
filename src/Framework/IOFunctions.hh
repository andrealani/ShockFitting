// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_IOFunctions_hh
#define ShockFitting_IOFunctions_hh

#ifdef SF_HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>

//----------------------------------------------------------------------------//


/// print an array
template <typename ARRAY>
static void print(std::ostream& out, const std::string& name, 
		  const ARRAY& v, int size, bool serial = false) 
{
  int rank = 0;
  
#ifdef SF_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (!serial) {
    out << "P" << rank << " => ";
  }
  
  if ((!serial) || rank == 0) {
    out << name << " = ";
    for (int w = 0; w < size; ++w) {out << v[w] << " ";}
    out << "\n";
  }
}

//----------------------------------------------------------------------------//

/// print a vector of std::pair
template <typename T1, typename T2> 
static void print(std::ostream& out, const std::string& name, 
		  const std::vector<std::pair<T1,T2> >& p, 
		  bool serial = false) 
{
  int rank = 0;
  
#ifdef SF_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (!serial) {
    out << "P" << rank << " => ";
  } 
  
  if ((!serial) || rank == 0) {
    out << name << " => ";
    for (unsigned i = 0; i < p.size(); ++i) {
      out << "[" << p[i].first << "," << p[i].second << "], ";
    }
    out << "\n";
  }
}

/// check integrity of the data
static void validate(bool ex, std::string errMsg, bool withRank = false) 
{
  int rank = 0;
#ifdef SF_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (!ex) {
    if (withRank) {std::cout << "on P" << rank << " ";}
    std::cout << errMsg << std::endl; abort();
  }
}

/// check if a file exists 
static bool fileExists(const char *filename)
{
  std::ifstream ifile(filename); return ifile;
}
#endif
