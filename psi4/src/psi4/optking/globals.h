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

/*! \file globals.h
    \ingroup optking
    \brief */

#ifndef _opt_globals_h_
#define _opt_globals_h_

const char* getIntcoFileName();
#define FILENAME_INTCO_DAT    getIntcoFileName()

#include <cstdio>

#include "package.h"
#if defined(OPTKING_PACKAGE_PSI)
  #include "psi4/psi4-dec.h"
#endif

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

namespace opt {
  EXTERN std::string psi_outfile;
  EXTERN FILE *qc_outfile;

  int read_natoms(void);
}


// symmetric matrix offset lookup array
#define IOFF_MAX 32641
namespace opt {
  EXTERN int *ioff;
}

// struct holding options/parameters for optking execution
#include "opt_params.h"
namespace opt {
  EXTERN OPT_PARAMS Opt_params;
}

// class for storage and manipulation of optimization step data
#include "opt_data.h"
namespace opt {
  EXTERN OPT_DATA *p_Opt_data;
}

#include "IRC_data.h"
namespace opt {
  EXTERN IRC_DATA *p_irc_data;
}

#include "opt_except.h"

#endif
