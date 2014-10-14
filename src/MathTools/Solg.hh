// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MathTools_Solg_hh
#define MathTools_Solg_hh

//--------------------------------------------------------------------------//

#include "MathTools/Array2D.hh"
#include <cstdlib>
#include <vector>

//--------------------------------------------------------------------------//

/// This class solves linear algebric system using Gauss Method.

template <class TYPE>
class Solg {
public:

  /// Constructor
  Solg() {}

  /// Destructor
  ~Solg() {}

  /// return vector of the unknows
  std::vector <TYPE> callSolg(Array2D <TYPE> A, std::vector <TYPE> B) {
    TYPE summ, pik, app;
    unsigned imax;
    std::vector <TYPE> X;
    X.resize(B.size());

    for (unsigned k=0; k<B.size()-1; k++) {
     imax = k;
     for (unsigned i1=k+1; i1<B.size(); i1++) {
      if (abs(A(i1,k)) > abs(A(imax,k))) {imax = i1;}
     }
     if (imax == k) {goto seven;}
     for (unsigned j1=k; j1<B.size(); j1++) {
      app = A(imax,j1);
      A(imax,j1) = A(k,j1);
      A(k,j1) = app;
     }
     app = B.at(k);
     B.at(k) = B.at(imax);
     B.at(imax) = app;
     seven:
          for (unsigned i=k+1; i<B.size(); i++) {
           pik = A(i,k)/A(k,k);
           A(i,k) = 0;
           B.at(i) = B.at(i)-pik*B.at(k);
           for (unsigned j=k+1; j<B.size(); j++) {
            A(i,j) = A(i,j) - pik * A(k,j);
           }
          }
    }

    X.at(B.size()-1) = B.at(B.size()-1)/A(B.size()-1,B.size()-1);
    for(int i=B.size()-2; i>=0; i--) {
     summ = 0;
     for (unsigned j=i+1; j<B.size(); j++) {
      summ = summ + A(i,j) * X.at(j);
     }
     X.at(i) = (B.at(i)-summ)/A(i,i);
    }
    return X;
 }

};

//--------------------------------------------------------------------------//

#endif  // MathTools_Solg_hh

