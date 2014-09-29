// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Log_hh
#define ShockFitting_Log_hh

#ifdef SF_HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>

//----------------------------------------------------------------------------//

/// enumerator holding all levels 
enum LogLevel {INFO = 0, VERBOSE = 1, DEBUG_MIN = 2, DEBUG_MAX = 3};

/// Singleton class holding the log level information
/// @author Andrea Lani
class LogSingleton {
public:
  
  /// get instance of this object
  static LogSingleton& getInstance() {static LogSingleton l; return l;}
  
  /// set the log level
  void setLevel(const LogLevel level) { m_level = level;}
  
  /// get the log level
  LogLevel getLevel() {return m_level;}
  
  /// get the log level name
  std::string getLevelName(const LogLevel level) {return m_levelToName.find(level)->second;}
  
  /// get the rank
  int getRank() 
  {
    int rank = 0;
#ifdef SF_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
  }
  
private:
  
  /// private constructor to make this class non-instantiable 
  LogSingleton() 
  {
    m_level = INFO;
    m_levelToName[INFO]      = "INFO";
    m_levelToName[VERBOSE]   = "VERBOSE";
    m_levelToName[DEBUG_MIN] = "DEBUG_MIN";
    m_levelToName[DEBUG_MAX] = "DEBUG_MAX";
  }
  
  /// private copy constructor to make this class non-instantiable 
  LogSingleton(const LogSingleton& l) {}
  
  /// log level
  LogLevel m_level;
  
  /// mapping between level and level name
  std::map<LogLevel, std::string> m_levelToName;
  
};

/// macro to print log info according to the selected level
#define LogToScreen(__level__, __content__)				\
  if ((unsigned)__level__ <= (unsigned) (LogSingleton::getInstance().getLevel())) \
    std::cout << "[P" << LogSingleton::getInstance().getRank() << "] [" \
	      << LogSingleton::getInstance().getLevelName((LogLevel)__level__) << "] " \
	      << __content__;

//----------------------------------------------------------------------------//

#endif
