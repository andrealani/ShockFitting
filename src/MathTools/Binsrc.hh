// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Binsrc_hh
#define Binsrc_hh

#include <vector>

//----------------------------------------------------------------------------//

/// This class finds number in sorted list (ascending)
/// with binary search
///

class Binsrc {
public:
 
  /// Constructor
  Binsrc();

  /// Constructor with kelem and klist assigned
  Binsrc (const unsigned, std::vector <int>);

  /// Destructor
  ~Binsrc();

  /// return ipos value
  int callBinsrc();

  /// return last value
  int getLast();

private: //data
 
  /// =-1: kelem not in table
  /// >=0: position in table
  int ipos;

  /// for ipos=-1, position behind which number belongs
  int last;

  /// number to be looked up
  int kelem;

  ///table vector in which the number is searched
  std::vector <int> klist;

};

//----------------------------------------------------------------------------//

#endif // Binsrc_hh

