/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_bin_psimrcc_debugging_h
#define _psi_src_bin_psimrcc_debugging_h
/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***********
****************************************************************/
#include "psi4/liboptions/liboptions.h"

namespace psi{ namespace psimrcc{

// #define NODEBUG
#ifdef NODEBUG
#define DEBUGGING(level,statements)
#else
#define DEBUGGING(level,statements) \
  if(debugging->is_level(level)){   \
    statements                      \
  }
#endif

#ifdef NODEBUG
#define START_TIMER(level,title)
#else
#define START_TIMER(level,title)             \
  Timer timer;                               \
  if(debugging->is_level(level)){            \
    outfile->Printf("\n  %-48s ...",title);  \
                             \
  }
#endif

#ifdef NODEBUG
#define END_TIMER(level)
#else
#define END_TIMER(level)                                    \
  if(debugging->is_level(level)){                           \
    outfile->Printf(" done. Timing %10.4f s",timer.get());  \
                                            \
  }
#endif

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class Debugging{
public:
  Debugging(Options &options);
  ~Debugging() {delete[] level;}

  bool is_level(int n) {return(n <= 10 ? level[n] : false);}
private:
  Options &options_;
  bool* level;
};

extern Debugging* debugging;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_debugging_h
