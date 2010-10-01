/*! \file globals.h
    \ingroup opt
    \brief */

#ifndef _opt_globals_h_
#define _opt_globals_h_

#include <cstdio>
#include "package.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

// get outfile pointer
#ifdef PSI4
#include <psi4-dec.h>
using psi::outfile;
#elif QCHEM4
namespace opt { EXTERN FILE *outfile; }
#endif

// set return type
#ifdef PSI4
typedef psi::PsiReturnType OptReturnType;
#define OptReturnEndloop (psi::Endloop)
#define OptReturnSuccess (psi::Success)
#elif QCHEM4
typedef int OptReturnType int;
#define OptReturnEndloop 1
#define OptReturnSuccess 0
#endif

// ioff array limit
#define IOFF_MAX 32641

// parameters for optimization and pointer for storage of optimization data
#include "opt_params.h"
#include "opt_data.h"
namespace opt {
  EXTERN OPT_PARAMS Opt_params;
  EXTERN OPT_DATA *p_Opt_data;

  EXTERN int *ioff;
}

#endif

