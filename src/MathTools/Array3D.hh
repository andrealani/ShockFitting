// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef Array3D_hh
#define Array3D_hh

//--------------------------------------------------------------------------//

#include <vector>
#include <assert.h>

//--------------------------------------------------------------------------//

/// This class defines an Array3D manipulationg one-dimensional objects
///

template <class TYPE>
class Array3D {
public:

  /// Construtor
  Array3D() {alreadyAlloc = false;}

  /// Constructor with a given size
  Array3D (const unsigned totsize) {size1 = totsize; alreadyAlloc=false;}

  /// Constructor with all sizes assigned
  Array3D(const unsigned dim1,
          const unsigned dim2,
          const unsigned dim3)
  {
    setSize1(dim1); setSize2(dim2);setSize3(dim3);
    ptr = new TYPE[size1*size2*size3]; ///create space for array
    for (unsigned i=0; i<size1*size2*size3; i++) { ptr[i] = 0;} ///initialize array
    alreadyAlloc = true;
  }

  ///Constructor with a given address
  Array3D (const unsigned dim1, const unsigned dim2, const unsigned dim3, TYPE* start)
  {
    setSize1(dim1); setSize2(dim2); setSize3(dim3);
    ptr = start;
    alreadyAlloc = false;
  }

  /// initialize array with all sizes assigned to 0
  void init(const unsigned dim1, const unsigned dim2, const unsigned dim3)
  {
    setSize1(dim1); setSize2(dim2);setSize3(dim3);
    ptr = new TYPE[size1*size2*size3]; //create space for array
    for (unsigned i=0; i<size1*size2*size3; i++) {ptr[i] = 0;}
  }

  /// initialize Array3D with a given address (sizes assigned)
  void address (const unsigned dim1, const unsigned dim2, const unsigned dim3, TYPE* start)
  {
    setSize1(dim1); setSize2(dim2); setSize3(dim3);
    ptr = start;
  }

  /// initialize Array3D with a given address (sizes not assigned)
  void address (TYPE* start) {ptr = start;}

  /// set array sizes
  void setSize1 (const unsigned dim) {size1=dim;}
  void setSize2 (const unsigned dim) {size2=dim;}
  void setSize3 (const unsigned dim) {size3=dim;}

  /// return array sizes
  unsigned getSize1() {return size1;}
  unsigned getSize2() {return size2;}
  unsigned getSize3() {return size3;}

  /// return total size
  unsigned size() {return size1*size2*size3;}

  /// overloading of the operator "()"
  TYPE & operator()(const unsigned i, const unsigned j, const unsigned k)
  {
    assert(i<size1); assert(j<size2); assert(k<size3);
    return ptr[k*size1*size2+j*size1+i];
  }

private: //data
  
  /// Verify allocation
  bool alreadyAlloc;

  /// Array3D sizes
  unsigned size1, size2, size3;

  /// pointer to first array element
  TYPE *ptr;
 
};

//--------------------------------------------------------------------------//

#endif //Array3D_hh

