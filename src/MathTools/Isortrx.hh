// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Isortrx_hh
#define Isortrx_hh

#include <vector>

//--------------------------------------------------------------------------//

/// This function performs an in-memory sort of the first N elements
/// of vector DATA, returning into vector INDEX the indices of elements
/// of DATA arranged in ascending order.  Thus,
///
///    DATA(INDEX(0)) will be the smallest number in vector DATA;
///    DATA(INDEX(N-1)) will be the largest number in DATA.
///
/// The original data is not physically rearranged.  The original order
/// of equal input values is not necessarily preserved.
///

class Isortrx {
public:

  /// Constructor
  Isortrx();

  /// Constructor
  Isortrx (std::vector <int>, const unsigned*);

  /// Destructor
  ~Isortrx();

  /// return vector INDEX with indices sorted in ascending order
  /// @ param INDEX : vector of DATA indices arranged in ascending order
  std::vector <int> callIsortrx();

private: //data

  /// dimension of vector DATA
  unsigned N;

  /// vector of index sorting
  std::vector <int> DATA;

};

//--------------------------------------------------------------------------//

#endif // Isortrx_hh
