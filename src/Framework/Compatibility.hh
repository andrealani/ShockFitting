// Copyright (C) 2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Compatibility_hh
#define ShockFitting_Compatibility_hh

// C++ to Fortran compatibility
#ifdef __STDC__
#define FSUBROUTINE(n_)  n_ ## _
#else
#define FSUBROUTINE(n_)  n_/**/_
#endif

//#ifdef CXX_NEEDS_FRIEND_TMPL_DECL
#  define LTGT <>
//#else
//#  define LTGT
//#endif

#endif
