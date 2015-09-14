// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/MovingAverageFilter.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

MovingAverageFilter::MovingAverageFilter()
{
}

//--------------------------------------------------------------------------//

MovingAverageFilter::~MovingAverageFilter()
{
}

//--------------------------------------------------------------------------//

void MovingAverageFilter::curveForwardSmoothing(vector<double> Y,
                                                int nbSmoothingPoints)
{
  int nbPoints = Y.size();

  Ynew.assign(nbPoints,0);

  // compute the new Y-ccordinates making an average of the Y-coordinates
  // of the closest points
  for(int IPOIN=0;IPOIN<nbPoints;IPOIN++) {
    
   if(IPOIN<=(nbPoints-nbSmoothingPoints)) {
    for(int ISMOOTHINGPOINT=0;
            ISMOOTHINGPOINT<nbSmoothingPoints;
            ISMOOTHINGPOINT++) {
     Ynew.at(IPOIN) += Y.at(IPOIN+ISMOOTHINGPOINT);
    }
   Ynew.at(IPOIN) /= nbSmoothingPoints;
   }
   else {
    for(int ISMOOTHINGPOINT=0;
            ISMOOTHINGPOINT<(nbPoints-IPOIN);
            ISMOOTHINGPOINT++) {
     Ynew.at(IPOIN) += Y.at(IPOIN+ISMOOTHINGPOINT);
    }
   Ynew.at(IPOIN) /= (nbPoints-IPOIN);
   }
  } 
}

//--------------------------------------------------------------------------//

