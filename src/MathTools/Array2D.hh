// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Array2D_hh
#define Array2D_hh

//--------------------------------------------------------------------------//

#include <vector>
#include <assert.h>

//--------------------------------------------------------------------------//

/// This class defines an Array2D manipulating one-dimensional objects
///

template <class TYPE>
class Array2D {
public:

 /// Default Constructor
  Array2D() {alreadyAlloc = false;}

 /// Constructor with a given size
 Array2D (const unsigned totsize) {nRows = totsize; alreadyAlloc = false;}

 /// Constructor with all assigned sizes
 Array2D(const unsigned nbRows, const unsigned nbColumns)
 {
 setnRows(nbRows); setnCols(nbColumns);
 ptr = new TYPE[nRows*nCols]; ///create space for the array
 for (unsigned i=0; i<nRows*nCols; i++) {ptr[i] = 0;}; ///initialize array
 alreadyAlloc = true;
 }

 /// Constructor
 /// set Array2D to a given address
 Array2D (const unsigned nbRows, const unsigned nbColumns, TYPE* start)
 {
  setnRows(nbRows); setnCols(nbColumns);
  ptr = start;
  alreadyAlloc = false;
 }

 /// Destructor
 ~Array2D(){if (alreadyAlloc) {delete ptr;}}

 /// Initialize array with given sizes to 0
 void init(const unsigned nbRows, const unsigned nbColumns)
 {
  setnRows(nbRows); setnCols(nbColumns);
  ptr = new TYPE[nRows*nCols];
  for (unsigned i=0; i<nRows*nCols; i++) {ptr[i] = 0;} //initialize array
  alreadyAlloc = true;
 }

 /// Set Array2D to an assigned address
 void address (const unsigned nbRows, const unsigned nbColumns, TYPE* start)
 {setnRows(nbRows); setnCols(nbColumns); ptr=start;}

 /// Set Array2D to an assigned address
 void address (TYPE* start) {ptr = start;}

 /// Set array sizes
 void setnRows (const unsigned dim) {nRows=dim;}
 void setnCols (const unsigned dim) {nCols=dim;}
 
 /// Return array sizes
 unsigned getnRows() const {return nRows;}
 unsigned getnCols() const {return nCols;}

 /// Return total size
 unsigned size() {return nRows*nCols;}
      
 /// Overloading of the operator()
 TYPE & operator()(const unsigned i, const unsigned j)
 {assert(i<nRows); assert(j<nCols); return ptr[j*nRows+i];}

 private:

  /// Verify allocation
  bool alreadyAlloc;

  /// Number of rows
  unsigned nRows;

  /// Number of columns
  unsigned nCols;

  /// Pointer to the first array element
  TYPE *ptr;

};

//--------------------------------------------------------------------------//

#endif //Array2D_hh

