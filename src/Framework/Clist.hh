// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Clist_hh
#define ShockFitting_Clist_hh

#include <vector>

#include "Framework/Compatibility.hh"
#include "SConfig/StringManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

template <class T> class Clist;
template <typename T> bool operator== (const T& value, const Clist<T>& c); 
template <typename T> bool operator!= (const T& value, const Clist<T>& c); 

// This class implements a list of values against which to compare
// @author Andrea Lani

template <typename T>
class Clist {
  // private definition of macro
#define L(a) m_v[a]=v##a
#define L1 L(0)
#define L2 L1;L(1)
#define L3 L2;L(2)
#define L4 L3;L(3)
#define L5 L4;L(4)
#define L6 L5;L(5)
#define L7 L6;L(6)
#define L8 L7;L(7)
#define L9 L8;L(8)
#define L10 L9;L(9)
#define L11 L10;L(10)
#define L12 L11;L(11)
#define L13 L12;L(12)
#define L14 L13;L(13)
#define L15 L14;L(14)
#define L16 L15;L(15)
#define L17 L16;L(16)
#define L18 L17;L(17)
#define L19 L18;L(18)
#define L20 L19;L(19)
#define L21 L20;L(20)
#define L22 L21;L(21)
#define L23 L22;L(22)
#define L24 L23;L(23)
#define L25 L24;L(24)
#define L26 L25;L(25)
#define L27 L26;L(26)
#define L28 L27;L(27)
#define L29 L28;L(28)
#define L30 L29;L(29)
#define L31 L30;L(30)
#define L32 L31;L(31)
#define L33 L32;L(32)
#define L34 L33;L(33)
#define L35 L34;L(34)
#define L36 L35;L(35)
#define L37 L36;L(36)
#define L38 L37;L(37)
#define L39 L38;L(38)
#define L40 L39;L(39)

#define A(a) T v##a
#define A1 A(0)
#define A2 A1,A(1)
#define A3 A2,A(2)
#define A4 A3,A(3)
#define A5 A4,A(4)
#define A6 A5,A(5)
#define A7 A6,A(6)
#define A8 A7,A(7)
#define A9 A8,A(8)
#define A10 A9,A(9)
#define A11 A10,A(10)
#define A12 A11,A(11)
#define A13 A12,A(12)
#define A14 A13,A(13)
#define A15 A14,A(14)
#define A16 A15,A(15)
#define A17 A16,A(16)
#define A18 A17,A(17)
#define A19 A18,A(18)
#define A20 A19,A(19)
#define A21 A20,A(20)
#define A22 A21,A(21)
#define A23 A22,A(22)
#define A24 A23,A(23)
#define A25 A24,A(24)
#define A26 A25,A(25)
#define A27 A26,A(26)
#define A28 A27,A(27)
#define A29 A28,A(28)
#define A30 A29,A(29)
#define A31 A30,A(30)
#define A32 A31,A(31)
#define A33 A32,A(32)
#define A34 A33,A(33)
#define A35 A34,A(34)
#define A36 A35,A(35)
#define A37 A36,A(36)
#define A38 A37,A(37)
#define A39 A38,A(38)
#define A40 A39,A(39)

#define MACLIST(a) Clist(A##a) {m_v.resize(a);L##a;}
  
public:

  friend bool operator== LTGT (const T& value, const Clist<T>& c); 
  friend bool operator!= LTGT (const T& value, const Clist<T>& c); 
  
  MACLIST(1)
  MACLIST(2)
  MACLIST(3)
  MACLIST(4)
  MACLIST(5)
  MACLIST(6)
  MACLIST(7)
  MACLIST(8)
  MACLIST(9)
  MACLIST(10)
  MACLIST(11)
  MACLIST(12)
  MACLIST(13)
  MACLIST(14)
  MACLIST(15)
  MACLIST(16)
  MACLIST(17)
  MACLIST(18)
  MACLIST(19)
  MACLIST(20)
  MACLIST(21)
  MACLIST(22)
  MACLIST(23)
  MACLIST(24)
  MACLIST(25)
  MACLIST(26)
  MACLIST(27)
  MACLIST(28)
  MACLIST(29)
  MACLIST(30)
  MACLIST(31)
  MACLIST(32)
  MACLIST(33)
  MACLIST(34)
  MACLIST(35)
  MACLIST(36)
  MACLIST(37)
  MACLIST(38)
  MACLIST(39)
  MACLIST(40)

private:
  // array storing the values
  std::vector<T> m_v;
};

//----------------------------------------------------------------------------//

template <typename T>
bool operator== (const T& value, const Clist<T>& c) 
{
  for (unsigned i = 0; i < c.m_v.size(); ++i) {
    if (c.m_v[i] == value) return true;
  }
  return false;
}

//----------------------------------------------------------------------------//

template <typename T>
bool operator!= (const T& value, const Clist<T>& c) 
{
  return !operator==(value,c);
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif
