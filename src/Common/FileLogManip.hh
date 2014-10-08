// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef FileLogManip_hh
#define FileLogManip_hh

//--------------------------------------------------------------------------//

#include <string>
#include <fstream>
#include <cstdlib>

//--------------------------------------------------------------------------//

/// This class defines FileLogManip, whose tasks are create a log file
/// with the name of the object which is referred to and memorize
/// input values in the log file using the overloading of operator ()

class FileLogManip {
public:

  /// Constructor
  FileLogManip() {}

  /// create file .log 
  void Open(std::string nameobj) { 
    namefile = nameobj;
    command = "touch " + namefile + ".log";
    system(command.c_str());
    command = "mv " + namefile + ".log ./log/";
    system(command.c_str());
    namefile = "./log/" + namefile + ".log";
    filelog.open(namefile.c_str());
  }
  /// close file .log
  void Close() {filelog.close();}

  /// Overloading of operator ()
  void operator() (std::string dummy) {
   filelog << dummy << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (T value) {
   filelog << value << " ";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (T value1, T value2) {
   filelog << value1 << ", " << value2 << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (const std::string dummy, T value) {
   filelog << dummy << " " << value << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (const std::string dummy, T value, std::string dummy2) {
   filelog << dummy << " " << value << " " << dummy2 << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (const std::string dummy, T value, 
                   std::string dummy2, std::string str) {
   filelog << dummy << " " << value << " " << dummy2 << " " << str << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (const std::string dummy1, T value1,
                    const std::string dummy2, T value2) {
   filelog << dummy1 << " " << value1 << " ";
   filelog << dummy2 << " " << value2 << "\n";}

  /// Overloading of operator ()
  template <typename T>
  void operator() (const std::string dummy1, T value1,
                    const std::string dummy2, T value2,
                    const std::string dummy3           ) {
   filelog << dummy1 << " " << value1 << " ";
   filelog << dummy2 << " " << value2 << " ";
   filelog << dummy3 <<  "\n";}

  /// Destructor
  ~FileLogManip() {}

private: //data

  /// name of created log file
  std::string namefile;

  /// dummy string used for "system" command
  std::string command;

  /// ofstream variable
  std::ofstream filelog;

};

#endif // FileLogManip


