// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <vector>
#include "Binsrc.hh"


#include<iostream>

//----------------------------------------------------------------------------//

Binsrc::Binsrc ()
{
}

//----------------------------------------------------------------------------//

Binsrc::~Binsrc()
{
}

//----------------------------------------------------------------------------//

Binsrc::Binsrc (const unsigned elem, std::vector <int> list)
{
  klist.resize(list.size());
  kelem = elem; klist = list; ipos=-1; last=0;
}

//----------------------------------------------------------------------------//

int Binsrc::callBinsrc ()
{
  int m,kpos;
  int nlist = klist.size();
  if ( kelem == klist.at(0) ){ipos = 0;}

  if ( nlist > 0) {
   kpos = 0;
  rep:
     m = (nlist+1)/2;
     last = kpos+m;
     if ( kelem < klist.at(last) ) {
      nlist = m;
      last = last-1;
      if ( (nlist) > 1) {goto rep;}
     }
     else if ( kelem > klist.at(last) ) {
      kpos = last;
      nlist = nlist-m;
      if (nlist>0) {goto rep;}
     }
     else {ipos = last;}
  }

  return ipos;
}

//----------------------------------------------------------------------------//

int Binsrc::getLast() {return last;}
