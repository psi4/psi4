
/*
 *  globals.h
 *  
 *
 *  Created by M.Saitow on 11/07/15.
 *  Copyright 2009 M.Saitow. All rights reserved.
 *
 */

#ifndef _psi_src_bin_adc_globals_h_
#define _psi_src_bin_adc_globals_h_

#include <ccfiles.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"

namespace psi { namespace adc {

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#define IOFF_MAX 32641 
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
	
//extern "C" {
//EXTERN FILE *infile, *outfile;
//EXTERN char *psi_file_prefix;
//}
EXTERN int *ioff;
EXTERN double **guess;
EXTERN double **MU_X;
EXTERN double **MU_Y;
EXTERN double **MU_Z;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
	
}}
#endif

