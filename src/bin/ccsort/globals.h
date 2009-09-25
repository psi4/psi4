/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>

namespace psi { namespace ccsort {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

extern "C" {
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
}

EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;

}} // namespace psi::ccsort
