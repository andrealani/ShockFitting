// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/Area.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

Area::Area()
{
}

//--------------------------------------------------------------------------//

Area::~Area()
{
}

//--------------------------------------------------------------------------//

void Area::computeArea(vector<double> X, vector <double> Y,
                       vector <unsigned> NODES)
{
  x = X; y = Y; nodes = NODES;
  x.resize(X.size());
  y.resize(Y.size());
  nodes.resize(NODES.size());

  unsigned NB = nodes.size();

  unsigned NNB; //logical copy of NB
  int ND;
  double x0, y0, DX1, DX2, DY1, DY2;
  double A; // partial sum of signed (and doubled) triangle areas

  NNB = NB;
  A = 0;
  if (NNB < 3) {goto two;}

  // initialization
  ND = nodes.at(0);
  x0 = x.at(ND);
  y0 = y.at(ND);
  ND = nodes.at(1);
  DX1 = x.at(ND)-x0;
  DY1 = y.at(ND)-y0;

  // loop on triangles (nodes[0],nodes[I],nodes[I+1]), I = 1, 3, ..., NB-2,
  // adding twice their signed areas to A
  for (unsigned I=2; I<NNB; I++) {
   ND = nodes.at(I);
   DX2 = x.at(ND) - x0;
   DY2 = y.at(ND) - y0;
   A = A + DX1*DY2 - DX2*DY1;
   DX1 = DX2;
   DY1 = DY2;
  }
  // A contains twice the signed area of the region
  two:
    area = A/2;

    return;
}

//--------------------------------------------------------------------------//
