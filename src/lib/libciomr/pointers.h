#ifndef _psi_src_lib_libciomr_pointers_h_
#define _psi_src_lib_libciomr_pointers_h_

#include "iomrparam.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

struct pointer {
     PSI_FPTR *wptr;
     };

EXTERN struct pointer ptr;
EXTERN int sector;
EXTERN time_t time_start, time_end;
#if !defined(SGI)
EXTERN struct tms total_tmstime;
#endif

#undef EXTERN

#ifdef __cplusplus
}
#endif

#endif /* header guard */
