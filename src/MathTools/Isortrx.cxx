// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include <vector>
#include "Isortrx.hh"

using namespace std;

//--------------------------------------------------------------------------//

Isortrx::Isortrx()
{
}

//--------------------------------------------------------------------------//

Isortrx::Isortrx (vector <int> data, const unsigned* n)
{
 DATA.resize(data.size());
 DATA = data;
 N = (*n);
}

//--------------------------------------------------------------------------//

Isortrx::~Isortrx()
{
}

//--------------------------------------------------------------------------//

vector <int> Isortrx::callIsortrx()
{
  int ISTK;
  unsigned L, R, I, J, P, INDEXP, INDEXT;
  int DATAP;
  vector <unsigned> LSTK;
  vector <unsigned> RSTK;
  vector <int> INDEX;

  INDEX.resize(N);

  LSTK.resize(31);
  RSTK.resize(31);


  const int M = 9;

  for (I=0; I<N; I++) { INDEX.at(I) = I;}

  if ( N < M) {goto nine;}

  ISTK = 0;
  L = 0;
  R = N-1;

  two:
    I = L;    J = R;

    P=(L+R)/2;
    INDEXP=INDEX.at(P);
    DATAP=DATA.at(INDEXP);

    if (DATA.at(INDEX.at(L)) > DATAP) {
     INDEX.at(P)=INDEX.at(L);
     INDEX.at(L)=INDEXP;
     INDEXP=INDEX.at(P);
     DATAP=DATA.at(INDEXP);
    }

    if (DATAP > DATA.at(INDEX.at(R))) {
     if (DATA.at(INDEX.at(L)) > DATA.at(INDEX.at(R))) {
      INDEX.at(P) = INDEX.at(L);
      INDEX.at(L) = INDEX.at(R);
     }
     else {
      INDEX.at(P) = INDEX.at(R);
     }
     INDEX.at(R) = INDEXP;
     INDEXP = INDEX.at(P);
     DATAP = DATA.at(INDEXP);
    }

   three:
        I = I+1;
        if (DATA.at(INDEX.at(I)) < DATAP) {goto three;}
    four:
        J = J-1;
        if (DATA.at(INDEX.at(J)) > DATAP) {goto four;}
     if (I < J) {
      INDEXT=INDEX.at(I);
      INDEX.at(I)=INDEX.at(J);
      INDEX.at(J)=INDEXT;
      goto three;
      }
     else {
      if ((R-J) >= (I-L) && (I-L) > M) {
       ISTK=ISTK+1;
       LSTK.insert(LSTK.begin()+ISTK, J+1);
       RSTK.insert(RSTK.begin()+ISTK, R);
       R=I-1;
      }
      else if ((I-L) > (R-J) && (R-J) > M) {
       ISTK=ISTK+1;
       LSTK.insert(LSTK.begin()+ISTK, L);
       RSTK.insert(RSTK.begin()+ISTK, I-1);
       L=J+1;
       }
      else if ((R-J) > M) {
       L=J+1;
      }
      else if ((I-L) > M) {
       R=I-1;
      }
      else {
       if (ISTK < 1) {goto nine;}
      L=LSTK.at(ISTK);
      R=RSTK.at(ISTK);
      ISTK=ISTK-1;
      }
      goto two;
     }

  nine:
    for (I=1; I< N; I++) {
     if (DATA.at(INDEX.at(I-1)) > DATA.at(INDEX.at(I))) {
      INDEXP = INDEX.at(I);
      DATAP=DATA.at(INDEXP);
      P=I-1;


  ninetwo:
       INDEX.at(P+1) = INDEX.at(P);
       P=P-1;
       if (P > 0) {
        if (DATA.at(INDEX.at(P)) > DATAP) {goto ninetwo;}
       }
       INDEX.at(P+1) = INDEXP;
      }
     }

    return INDEX;
}

//--------------------------------------------------------------------------//
