/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>

namespace psi { namespace cceom {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN void check_sum(char *term_lbl, int index, int irrep);

// #define EOM_DEBUG (0)
#define TIME_CCEOM (1)

#define H_IRR (0)
#define MAX(I,J) ((I>J) ? I : J)
#define MIN(I,J) ((I<J) ? I : J)

EXTERN FILE *outfile;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Eom_params eom_params;
EXTERN struct Local local;
EXTERN int ***dpd_dp;

}} // namespace psi::cceom
