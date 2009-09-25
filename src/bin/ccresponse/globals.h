/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>

namespace psi { namespace ccresponse {

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

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

}} // namespace psi::ccresponse
