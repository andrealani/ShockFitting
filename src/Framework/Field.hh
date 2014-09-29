// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Field_hh
#define ShockFitting_Field_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/Connectivity.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

enum {XX=0, YY=1, ZZ=2};  

/// This class defines a field. The actual field array MUST include 
/// both solution variables and spatial coordinates (x,y,z).  
/// 
/// @author Andrea Lani

class Field {
public:
  /// default constructor
  Field() : m_size(0), m_stride(0), m_conn(), m_values(NULL) {}
  
  /// constructor
  Field(const unsigned size, const unsigned stride, 
	const Connectivity& conn, double *const values) 
  {reset(size, stride, conn, values);}
  
  /// destructor
  ~Field() {}
  
  /// reset the object
  void reset(const unsigned size, const unsigned stride, 
	     const Connectivity& conn, double *const values) 
  {m_size = size; m_stride = stride; m_conn = conn; m_values = values;}
  
  /// get the size of the field
  unsigned getSize() const {return m_size;}
  
  /// get the stride (number of variables) of the field
  unsigned getStride() const {return m_stride;}
  
  /// get the element-state connectivity of the field
  const Connectivity& getConnectivity() const {return m_conn;}
  
  /// get the field array
  double* getArray() const {return m_values;}
  
private:
  
  /// number of degrees of freedom
  unsigned m_size;
  
  /// stride for the degrees of freedom (e.g. number of equations)
  unsigned m_stride;
  
  /// connectivity 
  Connectivity m_conn;
  
  /// array of values
  double* m_values;
  
};
  
//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_Field_hh
