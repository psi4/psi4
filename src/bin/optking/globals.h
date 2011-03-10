/*! \file globals.h
    \ingroup optking
    \brief */

#ifndef _opt_globals_h_
#define _opt_globals_h_

#define FILENAME_INTCO_DAT    "intco.dat"
#define FILENAME_OUTPUT_DAT   "output.dat"     // not used by PSI
#define FILENAME_CARTESIAN_H  "psi.file15.dat" // used by PSI
#define FILENAME_GEOM_GRAD_IN "psi.file11.dat" // used by PSI

#include <cstdio>
#include "package.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#if defined(OPTKING_PACKAGE_PSI)
  #include <psi4-dec.h>
#endif

namespace opt {
  EXTERN FILE *outfile;

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

#include "opt_except.h"

#endif

