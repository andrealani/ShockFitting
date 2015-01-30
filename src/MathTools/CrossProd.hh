// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef CrossProd_hh
#define CrossProd_hh

//--------------------------------------------------------------------------//

/// This class defines a CrossProd, whose task is to compute
/// the cross product between two vectors

template <class TYPE>
class CrossProd {
public:

  /// Constructor
  CrossProd() {}

  /// Desctructor
  ~CrossProd() {}

  void computeCrossProd(std::vector <TYPE> a, std::vector <TYPE> b)
  {
    if(a.size()!=b.size()) {
     std::cout << "CrossProd::error in dimension size\n";
    }
    c.resize(a.size());
    c.at(0) = a.at(1) * b.at(2) - a.at(2) * b.at(1);
    c.at(1) = a.at(2) * b.at(0) - a.at(0) * b.at(2);
    c.at(2) = a.at(0) * b.at(1) - a.at(1) * b.at(0);
  }

  std::vector <TYPE> getCrossProd() const { return c; }

private:

  /// vector of the cross product
  std::vector <TYPE> c;
};

//--------------------------------------------------------------------------//

#endif // CrossProd_hh 
