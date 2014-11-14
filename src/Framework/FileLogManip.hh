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
    system("if [ ! -d log ]; then  mkdir log ; fi");
    namefile = nameobj;
    command = "touch " + namefile + ".log";
    system(command.c_str());
    command = "mv " + namefile + ".log ./log/";
    system(command.c_str());
    namefile = "./log/" + namefile + ".log";
    filelog.open(namefile.c_str());
    filelog.precision(16);
  }
  /// close file .log
  void Close() {filelog.close();}

  /// Overloading of operator ()
  template <typename T>
  void operator() (T dummy) {
   filelog << dummy << "\n";}

  /// Overloading of operator ()
  template <typename T1, typename T2>
  void operator() (T1 dummy1, T2 dummy2) {
   filelog << dummy1 << " " << dummy2;}

  /// Overloading of operator ()
  template <typename T1, typename T2, typename T3>
  void operator() (T1 dummy1, T2 dummy2, T3 dummy3) {
   filelog << dummy1 << " " << dummy2 << " " << dummy3;}

  /// Overloading of operator ()
  template <typename T1, typename T2, typename T3, typename T4>
  void operator() (T1 dummy1, T2 dummy2, T3 dummy3, T4 dummy4) {
   filelog << dummy1 << " " << dummy2 << " ";
   filelog << dummy3 << " " << dummy4;}

  /// Overloading of operator ()
  template <typename T1, typename T2, typename T3, 
              typename T4, typename T5>
  void operator() (T1 dummy1, T2 dummy2, T3 dummy3,
                   T4 dummy4, T5 dummy5) {
   filelog << dummy1 << " " << dummy2 << " " << dummy3;
   filelog << " " << dummy4 << " " << dummy5;}

  /// Overloading of operator ()
  template <typename T1, typename T2, typename T3,
            typename T4,typename T5, typename T6>
  void operator () (T1 dummy1, T2 dummy2, T3 dummy3,
                    T4 dummy4, T5 dummy5, T6 dummy6) {
   filelog << dummy1 << " " << dummy2 << " ";
   filelog << dummy3 << " " << dummy4 << " ";
   filelog << dummy5 << " " << dummy6;}


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


