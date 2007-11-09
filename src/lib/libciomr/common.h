#ifndef _psi_src_lib_libciomr_common_h_
#define _psi_src_lib_libciomr_common_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN FILE *intfile, *outfile;

#ifdef __cplusplus
}
#endif

#endif /* header guard */
