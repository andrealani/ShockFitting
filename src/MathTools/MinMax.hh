// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MinMax_hh
#define MinMax_hh

//--------------------------------------------------------------------------//

/// This class defines a MinMax, whose task is to find the minimum or
/// the maximum value bewteen several elements

template <class TYPE>
class MinMax {
public:

  /// Costructor
  MinMax() {}

  /// Destructor
  ~MinMax() {}


  /// return minimum value
  TYPE min(TYPE value1, TYPE value2, TYPE value3) {
   minvalue = value1;
   if (value2<minvalue) { minvalue = value2; }
   if (value3<minvalue) { minvalue = value3; }
   return minvalue;
  }

  /// return maximum value
  TYPE max(TYPE value1, TYPE value2, TYPE value3) {
   maxvalue = value1;
   if (value2>maxvalue) { maxvalue = value2; }
   if (value3>maxvalue) { maxvalue = value3; }
   return maxvalue;
  }

private: // data

  TYPE minvalue;

  TYPE maxvalue;
};

//--------------------------------------------------------------------------//

#endif // MinMax_hh
